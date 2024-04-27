use axiom_eth::rlc::{
    chip::RlcChip,
    circuit::{builder::RlcCircuitBuilder, instructions::RlcCircuitInstructions, RlcCircuitParams},
};
use halo2_base::{
    gates::{circuit::BaseCircuitParams, GateInstructions, RangeChip, RangeInstructions},
    utils::ScalarField,
    QuantumCell::Constant,
};
use serde::Deserialize;

use crate::{
    constants::{
        sk_enc_constants_1024_1x27_65537::S_BOUND,
        sk_enc_constants_4096_2x55_65537::{E_BOUND, K0IS, K1_BOUND, N, QIS, R1_BOUNDS, R2_BOUNDS},
    },
    poly::{Poly, PolyAssigned},
};

/// Helper function to define the parameters of the RlcCircuit. This is a non-optimized configuration that makes use of a single advice column. Use this for testing purposes only.
pub fn test_params() -> RlcCircuitParams {
    RlcCircuitParams {
        base: BaseCircuitParams {
            k: 21,
            num_advice_per_phase: vec![1, 1],
            num_fixed: 1,
            num_lookup_advice_per_phase: vec![0, 1],
            lookup_bits: Some(8),
            num_instance_columns: 0,
        },
        num_rlc_columns: 1,
    }
}
/// `BfvSkEncryptionCircuit` is a circuit that checks the correct formation of a ciphertext resulting from BFV secret key encryption
/// All the polynomials coefficients and scalars are normalized to be in the range `[0, p)` where p is the modulus of the prime field of the circuit
///
/// pk_q1 = ( pk0_qi , pk1_qi )=( [ai*s + E] , -ai )
/// # Parameters:
/// * `pko_qi`: publicly polynomial created by secret polynomial ([ai*s + E] )
/// * `u`: secret polynomial, sampled from ternary distribution.
/// * `e0`: error polynomial, sampled from discrete Gaussian distribution.
/// * `k1`: scaled message polynomial.
/// * `r2is`: list of r2i polynomials for each i-th CRT basis .
/// * `r1is`: list of r1i polynomials for each CRT i-th CRT basis.
/// * `ct0is`: list of ct0i (first component of the ciphertext cti) polynomials for each CRT i-th CRT basis.

#[derive(Deserialize, Clone)]
pub struct BfvPkEncryptionCircuit {
    pk0_qi: Vec<String>,
    u: Vec<String>,
    e0: Vec<String>,
    k1: Vec<String>,
    r2is: Vec<Vec<String>>,
    r1is: Vec<Vec<String>>,
    ct0is: Vec<Vec<String>>,
}

/// Payload returned by the first phase of the circuit to be reused in the second phase
pub struct Payload<F: ScalarField> {
    pk0_qi_assigned: PolyAssigned<F>,
    u_assigned: PolyAssigned<F>,
    e0_assigned: PolyAssigned<F>,
    k1_assigned: PolyAssigned<F>,
    r2is_assigned: Vec<PolyAssigned<F>>,
    r1is_assigned: Vec<PolyAssigned<F>>,
    ct0is: Vec<Vec<String>>,
}

impl<F: ScalarField> RlcCircuitInstructions<F> for BfvPkEncryptionCircuit {
    type FirstPhasePayload = Payload<F>;

    /// #### Phase 0
    /// TODO: Figure out whether to assign the pk0_qi in phase0 or not , verifieer knows the pk(doubt)

    /// In this phase, the polynomials for each matrix $S_i$ are assigned to the circuit. Namely:
    /// * polynomials `u`, `e0`, `k1`, `pk0_qi` are assigned to the witness table. This has to be done only once as these polynomial are common to each $S_i$ matrix
    /// * polynomials `r1i`, `r2i` are assigned to the witness table for each $S_i$ matrix

    /// Witness values are element of the finite field $\mod{p}$. Negative coefficients $-z$ are assigned as field elements $p - z$.

    /// At the end of phase 0, the witness generated so far is interpolated into a polynomial and committed by the prover. The hash of this commitment is used as challenge and will be used as a source of randomness $\gamma$ in Phase 1. This feature is made available by Halo2 [Challenge API](https://hackmd.io/@axiom/SJw3p-qX3).

    fn virtual_assign_phase0(
        &self,
        builder: &mut RlcCircuitBuilder<F>,
        _: &RangeChip<F>,
    ) -> Self::FirstPhasePayload {
        let ctx = builder.base.main(0);

        let pk0_qi = Poly::<F>::new(self.pk0_qi.clone());
        let pk0_qi_assigned = PolyAssigned::new(ctx, pk0_qi);

        let u = Poly::<F>::new(self.u.clone());
        let u_assigned = PolyAssigned::new(ctx, u);

        let e0 = Poly::<F>::new(self.e0.clone());
        let e0_assigned = PolyAssigned::new(ctx, e0);

        let k1 = Poly::<F>::new(self.k1.clone());
        let k1_assigned = PolyAssigned::new(ctx, k1);

        let r1is_assigned = self
            .r1is
            .iter()
            .map(|r1is| {
                let r1is = Poly::<F>::new(r1is.clone());
                PolyAssigned::new(ctx, r1is)
            })
            .collect::<Vec<_>>();

        let r2is_assigned = self
            .r2is
            .iter()
            .map(|r2is| {
                let r2is = Poly::<F>::new(r2is.clone());
                PolyAssigned::new(ctx, r2is)
            })
            .collect::<Vec<_>>();

        Payload {
            pk0_qi_assigned,
            u_assigned,
            e0_assigned,
            k1_assigned,
            r2is_assigned,
            r1is_assigned,
            ct0is: self.ct0is.clone(),
        }
    }

    /// #### Phase 1

    /// In this phase, the following two core constraints are enforced:

    /// - The coefficients of $S_i$ are in the expected range.  
    /// - $P_i(\gamma) \times S_i(\gamma) =Ct_{0,i}(\gamma)$

    /// ##### Assignment
    /// * Assign evaluations to the circuit:  `ct0i(gamma)` for each $P_i$ matrix
    /// * Assign `cyclo(gamma)` to the circuit. This has to be done only once as the cyclotomic polynomial is common to each $U_i$ matrix
    /// * Expose `ct0i(gamma)` for each $U_i$ matrix
    /// * Expose `cyclo(gamma)` as public input

    /// NOT SURE ABOUT THIS , BUT AT PRESENT CONSIDERISING IT AS UNKNOWN TO VERIFIER
    /// Since the polynomials `cyclo`, `pko_qi` and `ct0i` are known to the verifier, the evaluation at $\gamma$ doesn't need to be constrained inside the circuit. Instead, this can be safely be performed (and verified) outside the circuit.

    /// ##### Range Check

    /// The coefficients of the private polynomials from each $i$-th matrix $S_i$ are checked to be in the correct range
    /// * Range check polynomials `u`, `e0`, `k1`. This has to be done only once as these polynomial are common to each $S_i$ matrix
    /// * Range check polynomials `r1i`, `r2i` for each $S_i$ matrix

    /// Since negative coefficients `-z` are assigned as `p - z` to the circuit, this might result in very large coefficients. Performing the range check on such large coefficients requires large lookup tables. To avoid this, the coefficients (both negative and positive) are shifted by a constant to make them positive and then perform the range check.

    /// ##### Evaluation at $\gamma$ Constraint

    /// Contrary to the polynomials `cyclo`, `pk0_qi` and `ct0i`, the polynomials belonging to each $S_i$ matrix  are not known by the verifier. Therefore, their evaluation at $\gamma$ must be constrained inside the circuit.

    /// * Constrain the evaluation of the polynomials `u`, `e0`, `k1` at $\gamma$. This has to be done only once as these polynomial are common to each $S_i$ matrix
    /// * Constrain the evaluation of the polynomials `r1i`, `r2i` at $\gamma$ for each $S_i$ matrix

    /// ##### Correct Encryption Constraint

    /// It is needed to prove that $P_i(\gamma) \times S_i(\gamma) =Ct_{0,i}(\gamma)$. This can be rewritten as `ct0i = ct0i_hat + r1i * qi + r2i * cyclo`, where `ct0i_hat = pk0_qi * u + e0*1 + k1 * k0i`.

    /// This constrained is enforced by proving that `LHS(gamma) = RHS(gamma)`. According to the Schwartz-Zippel lemma, if this relation between polynomial when evaluated at a random point holds true, then then the polynomials are identical with high probability. Note that `qi` and `k0i` (for each $U_i$ matrix) are constants to the circuit encoded during key generation.
    /// * Constrain that `ct0i(gamma) = ai(gamma) * s(gamma) + e(gamma) + k1(gamma) * k0i + r1i(gamma) * qi + r2i(gamma) * cyclo(gamma)` for each $i$-th CRT basis
    ///

    fn virtual_assign_phase1(
        builder: &mut RlcCircuitBuilder<F>,
        range: &RangeChip<F>,
        rlc: &RlcChip<F>,
        payload: Self::FirstPhasePayload,
    ) {
        let Payload {
            pk0_qi_assigned,
            u_assigned,
            e0_assigned,
            k1_assigned,
            r2is_assigned,
            r1is_assigned,
            ct0is,
        } = payload;

        // ASSIGNMENT

        let (ctx_gate, ctx_rlc) = builder.rlc_ctx_pair();
        let gamma = *rlc.gamma();

        let mut ct0is_at_gamma_assigned = vec![];
        let mut qi_constants = vec![];
        let mut k0i_constants = vec![];

        for z in 0..ct0is.len() {
            let ct0i = Poly::<F>::new(ct0is[z].clone());
            let ct0i_at_gamma = ct0i.eval(gamma);
            let ct0i_at_gamma_assigned = ctx_gate.load_witness(ct0i_at_gamma);
            ct0is_at_gamma_assigned.push(ct0i_at_gamma_assigned);

            let qi_constant = Constant(F::from_str_vartime(QIS[z]).unwrap());
            qi_constants.push(qi_constant);

            let k0i_constant = Constant(F::from_str_vartime(K0IS[z]).unwrap());
            k0i_constants.push(k0i_constant);
        }

        // cyclo poly is equal to x^N + 1
        let cyclo_at_gamma = gamma.pow_vartime([N as u64]) + F::from(1);
        let cyclo_at_gamma_assigned = ctx_gate.load_witness(cyclo_at_gamma);

        u_assigned.range_check(ctx_gate, range, S_BOUND);
        e0_assigned.range_check(ctx_gate, range, E_BOUND);
        k1_assigned.range_check(ctx_gate, range, K1_BOUND);

        for z in 0..ct0is.len() {
            r2is_assigned[z].range_check(ctx_gate, range, R2_BOUNDS[z]);
            r1is_assigned[z].range_check(ctx_gate, range, R1_BOUNDS[z]);
        }

        let u_at_gamma = u_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);
        let e0_at_gamma = e0_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);
        let k1_at_gamma = k1_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);
        let pk0_qi_at_gamma = pk0_qi_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);

        let gate = range.gate();

        // For each `i` Prove that LHS(gamma) = RHS(gamma)
        // pk_u = pk0_qi(gamma) * u(gamma) + e0(gamma)
        // LHS = ct0i(gamma)
        // RHS = pk_u  + k1(gamma) * k0i + r1i(gamma) * qi + r2i(gamma) * cyclo(gamma)

        let pk_u = gate.mul_add(ctx_gate, pk0_qi_at_gamma, u_at_gamma, e0_at_gamma);

        for z in 0..ct0is.len() {
            let r1i_at_gamma = r1is_assigned[z].enforce_eval_at_gamma(ctx_rlc, rlc);
            let r2i_at_gamma = r2is_assigned[z].enforce_eval_at_gamma(ctx_rlc, rlc);

            // CORRECT ENCRYPTION CONSTRAINT

            // rhs = rhs + k1(gamma) * k0i
            let rhs = gate.mul_add(ctx_gate, k1_at_gamma, k0i_constants[z], pk_u);

            // rhs = rhs + r1i(gamma) * qi
            let rhs = gate.mul_add(ctx_gate, r1i_at_gamma, qi_constants[z], rhs);

            // rhs = rhs + r2i(gamma) * cyclo(gamma)
            let rhs = gate.mul_add(ctx_gate, r2i_at_gamma, cyclo_at_gamma_assigned, rhs);
            let lhs = ct0is_at_gamma_assigned[z];

            // LHS(gamma) = RHS(gamma)
            let res = gate.is_equal(ctx_gate, lhs, rhs);
            gate.assert_is_const(ctx_gate, &res, &F::from(1));
        }
    }
}
