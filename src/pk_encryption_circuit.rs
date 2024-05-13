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
    constants::pk_enc_constants::pk_enc_constants_1024_15x60_65537::{
        E_BOUND, K0IS, K1_BOUND, N, P1_BOUNDS, P2_BOUNDS, PK_BOUND, QIS, R1_BOUNDS, R2_BOUNDS,
        U_BOUND,
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
/// * `pk0_qi`: publicly polynomial created by secret polynomial ([ai*s + E] )
/// * `pk1_qi`: publicly polynomial created by polynomial (-[ai])
/// * `u`: secret polynomial, sampled from ternary distribution.
/// * `e0`: error polynomial, sampled from discrete Gaussian distribution.
/// * `e1`: error polynomial, sampled from discrete Gaussian distribution.
/// * `k1`: scaled message polynomial.
/// * `r2is`: list of r2i polynomials for each i-th CRT basis .
/// * `r1is`: list of r1i polynomials for each CRT i-th CRT basis.
/// * `p2is`: list of p2i polynomials for each i-th CRT basis .
/// * `ct0is`: list of ct0i (first component of the ciphertext cti) polynomials for each CRT i-th CRT basis.

#[derive(Deserialize, Clone)]
pub struct BfvPkEncryptionCircuit {
    pk0_qi: Vec<Vec<String>>,
    pk1_qi: Vec<Vec<String>>,
    u: Vec<String>,
    e0: Vec<String>,
    e1: Vec<String>,
    k1: Vec<String>,
    r2is: Vec<Vec<String>>,
    r1is: Vec<Vec<String>>,
    p2is: Vec<Vec<String>>,
    p1is: Vec<Vec<String>>,
    ct0is: Vec<Vec<String>>,
    ct1is: Vec<Vec<String>>,
}

/// Payload returned by the first phase of the circuit to be reused in the second phase
pub struct Payload<F: ScalarField> {
    pk0_qi_assigned: Vec<PolyAssigned<F>>,
    pk1_qi_assigned: Vec<PolyAssigned<F>>,
    u_assigned: PolyAssigned<F>,
    e0_assigned: PolyAssigned<F>,
    e1_assigned: PolyAssigned<F>,
    k1_assigned: PolyAssigned<F>,
    r2is_assigned: Vec<PolyAssigned<F>>,
    r1is_assigned: Vec<PolyAssigned<F>>,
    p2is_assigned: Vec<PolyAssigned<F>>,
    p1is_assigned: Vec<PolyAssigned<F>>,
    ct0is: Vec<Vec<String>>,
    ct1is: Vec<Vec<String>>,
}

impl<F: ScalarField> RlcCircuitInstructions<F> for BfvPkEncryptionCircuit {
    type FirstPhasePayload = Payload<F>;

    /// #### Phase 0

    /// In this phase, the polynomials for each matrix $S_i$ are assigned to the circuit. Namely:
    /// * polynomials `u`,'e1, `e0`, `k1`, `pk0_qi`,`pk1_q1` are assigned to the witness table. This has to be done only once as these polynomial are common to each $S_i$ matrix
    /// * polynomials `r1i`, `r2i`,`p1i`,`p2i` are assigned to the witness table for each $S_i$ matrix

    /// Witness values are element of the finite field $\mod{p}$. Negative coefficients $-z$ are assigned as field elements $p - z$.

    /// At the end of phase 0, the witness generated so far is interpolated into a polynomial and committed by the prover. The hash of this commitment is used as challenge and will be used as a source of randomness $\gamma$ in Phase 1. This feature is made available by Halo2 [Challenge API](https://hackmd.io/@axiom/SJw3p-qX3).

    fn virtual_assign_phase0(
        &self,
        builder: &mut RlcCircuitBuilder<F>,
        _: &RangeChip<F>,
    ) -> Self::FirstPhasePayload {
        let ctx = builder.base.main(0);

        let pk0_qi = self
            .pk0_qi
            .iter()
            .map(|pk0_qi| Poly::<F>::new(pk0_qi.clone()))
            .collect::<Vec<_>>();
        let pk0_qi_assigned = pk0_qi
            .into_iter()
            .map(|pk0_qi| PolyAssigned::new(ctx, pk0_qi))
            .collect::<Vec<_>>();

        let pk1_qi = self
            .pk1_qi
            .iter()
            .map(|pk1_qi| Poly::<F>::new(pk1_qi.clone()))
            .collect::<Vec<_>>();
        let pk1_qi_assigned = pk1_qi
            .into_iter()
            .map(|pk1_qi| PolyAssigned::new(ctx, pk1_qi))
            .collect::<Vec<_>>();

        let u = Poly::<F>::new(self.u.clone());
        let u_assigned = PolyAssigned::new(ctx, u);

        let e0 = Poly::<F>::new(self.e0.clone());
        let e0_assigned = PolyAssigned::new(ctx, e0);

        let e1 = Poly::<F>::new(self.e1.clone());
        let e1_assigned = PolyAssigned::new(ctx, e1);

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

        let p1is_assigned = self
            .p1is
            .iter()
            .map(|p1is| {
                let p1is = Poly::<F>::new(p1is.clone());
                PolyAssigned::new(ctx, p1is)
            })
            .collect::<Vec<_>>();

        let p2is_assigned = self
            .p2is
            .iter()
            .map(|p2is| {
                let p2is = Poly::<F>::new(p2is.clone());
                PolyAssigned::new(ctx, p2is)
            })
            .collect::<Vec<_>>();

        Payload {
            pk0_qi_assigned,
            pk1_qi_assigned,
            u_assigned,
            e0_assigned,
            e1_assigned,
            k1_assigned,
            r2is_assigned,
            r1is_assigned,
            p2is_assigned,
            p1is_assigned,
            ct0is: self.ct0is.clone(),
            ct1is: self.ct1is.clone(),
        }
    }

    /// #### Phase 1

    /// In this phase, the following two core constraints are enforced:

    /// - The coefficients of $S_i$ are in the expected range.  
    /// - $P_i(\gamma) \times S_i(\gamma) =Ct_{0,i}(\gamma)$

    /// ##### Assignment
    /// * Assign evaluations to the circuit:  `ct0i(gamma)` and `ct1i(gamma)` for each $P_i$ matrix
    /// * Assign `cyclo(gamma)` to the circuit. This has to be done only once as the cyclotomic polynomial is common to each $U_i$ matrix
    /// * Expose `ct0i(gamma)` and `ct1i(gamma)` for each $U_i$ matrix
    /// * Expose `cyclo(gamma)` as public input

    /// Since the polynomials `cyclo`, `pk0_qi` ,`pk1_qi` ,`ct1i` and `ct0i` are known to the verifier, the evaluation at $\gamma$ doesn't need to be constrained inside the circuit. Instead, this can be safely be performed (and verified) outside the circuit.

    /// ##### Range Check

    /// The coefficients of the private polynomials from each $i$-th matrix $S_i$ are checked to be in the correct range
    /// * Range check polynomials `u`, `e0`,`e1`,`k1`. This has to be done only once as these polynomial are common to each $S_i$ matrix
    /// * Range check polynomials `r1i`, `r2i` for each $S_i$ matrix
    /// * Range check polynomials `p1i`, `p2i` for each $S_i$ matrix

    /// Since negative coefficients `-z` are assigned as `p - z` to the circuit, this might result in very large coefficients. Performing the range check on such large coefficients requires large lookup tables. To avoid this, the coefficients (both negative and positive) are shifted by a constant to make them positive and then perform the range check.

    /// ##### Evaluation at $\gamma$ Constraint

    /// Contrary to the polynomials `cyclo`, `pk0_qi`,`pk1_qi` and `ct0i`, the polynomials belonging to each $S_i$ matrix  are not known by the verifier. Therefore, their evaluation at $\gamma$ must be constrained inside the circuit.

    /// * Constrain the evaluation of the polynomials `u`, `e0`, `e1`, `k1` at $\gamma$. This has to be done only once as these polynomial are common to each $S_i$ matrix
    /// * Constrain the evaluation of the polynomials `r1i`, `r2i` at $\gamma$ for each $S_i$ matrix
    /// * Constrain the evaluation of the polynomials `p1i`, `p2i` at $\gamma$ for each $S_i$ matrix

    /// ##### Correct Encryption Constraint

    /// It is needed to prove that $P_i(\gamma) \times S_i(\gamma) =Ct_{0,i}(\gamma)$. This can be rewritten as `ct0i = ct0i_hat + r1i * qi + r2i * cyclo`, where `ct0i_hat = pk0_qi * u + e0 + k1 * k0i`.

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
            pk1_qi_assigned,
            u_assigned,
            e0_assigned,
            e1_assigned,
            k1_assigned,
            r2is_assigned,
            r1is_assigned,
            p2is_assigned,
            p1is_assigned,
            ct0is,
            ct1is,
        } = payload;

        // ASSIGNMENT

        let (ctx_gate, ctx_rlc) = builder.rlc_ctx_pair();
        let gamma = *rlc.gamma();

        let mut ct0is_at_gamma_assigned = vec![];
        let mut ct1is_at_gamma_assigned = vec![];
        let mut qi_constants = vec![];
        let mut k0i_constants = vec![];

        for z in 0..ct0is.len() {
            let ct0i = Poly::<F>::new(ct0is[z].clone());
            let ct0i_at_gamma = ct0i.eval(gamma);
            let ct0i_at_gamma_assigned = ctx_gate.load_witness(ct0i_at_gamma);
            ct0is_at_gamma_assigned.push(ct0i_at_gamma_assigned);

            let ct1i = Poly::<F>::new(ct1is[z].clone());
            let ct1i_at_gamma = ct1i.eval(gamma);
            let ct1i_at_gamma_assigned = ctx_gate.load_witness(ct1i_at_gamma);
            ct1is_at_gamma_assigned.push(ct1i_at_gamma_assigned);

            let qi_constant = Constant(F::from_str_vartime(QIS[z]).unwrap());
            qi_constants.push(qi_constant);

            let k0i_constant = Constant(F::from_str_vartime(K0IS[z]).unwrap());
            k0i_constants.push(k0i_constant);
        }

        // cyclo poly is equal to x^N + 1
        let cyclo_at_gamma = gamma.pow_vartime([N as u64]) + F::from(1);
        let cyclo_at_gamma_assigned = ctx_gate.load_witness(cyclo_at_gamma);

        u_assigned.range_check(ctx_gate, range, U_BOUND);
        e0_assigned.range_check(ctx_gate, range, E_BOUND);
        e1_assigned.range_check(ctx_gate, range, E_BOUND);
        k1_assigned.range_check(ctx_gate, range, K1_BOUND);

        let _ = pk0_qi_assigned
            .iter()
            .enumerate()
            .map(|(i, pk_assigned)| pk_assigned.range_check(ctx_gate, range, PK_BOUND[i]));

        let _ = pk1_qi_assigned
            .iter()
            .enumerate()
            .map(|(i, pk_assigned)| pk_assigned.range_check(ctx_gate, range, PK_BOUND[i]));

        for z in 0..ct0is.len() {
            r2is_assigned[z].range_check(ctx_gate, range, R2_BOUNDS[z]);
            r1is_assigned[z].range_check(ctx_gate, range, R1_BOUNDS[z]);
            p2is_assigned[z].range_check(ctx_gate, range, P2_BOUNDS[z]);
            p1is_assigned[z].range_check(ctx_gate, range, P1_BOUNDS[z]);
        }

        let u_at_gamma = u_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);
        let e0_at_gamma = e0_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);
        let e1_at_gamma = e1_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);
        let k1_at_gamma = k1_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);
        let pk0_qi_at_gamma = pk0_qi_assigned
            .iter()
            .map(|pk_assigned| pk_assigned.enforce_eval_at_gamma(ctx_rlc, rlc))
            .collect::<Vec<_>>();
        let pk1_qi_at_gamma = pk1_qi_assigned
            .iter()
            .map(|pk_assigned| pk_assigned.enforce_eval_at_gamma(ctx_rlc, rlc))
            .collect::<Vec<_>>();
        let gate = range.gate();

        // For each `i` Prove that LHS(gamma) = RHS(gamma)
        // pk0_u = pk0_qi(gamma) * u(gamma) + e0(gamma)
        // LHS = ct0i(gamma)
        // RHS = pk0_u  + k1(gamma) * k0i + r1i(gamma) * qi + r2i(gamma) * cyclo(gamma)

        for z in 0..ct0is.len() {
            let pk0_u = gate.mul_add(ctx_gate, pk0_qi_at_gamma[z], u_at_gamma, e0_at_gamma);
            let r1i_at_gamma = r1is_assigned[z].enforce_eval_at_gamma(ctx_rlc, rlc);
            let r2i_at_gamma = r2is_assigned[z].enforce_eval_at_gamma(ctx_rlc, rlc);

            // CORRECT ENCRYPTION CONSTRAINT

            // rhs = pk0_u + k1(gamma) * k0i
            let rhs = gate.mul_add(ctx_gate, k1_at_gamma, k0i_constants[z], pk0_u);

            // rhs = rhs + r1i(gamma) * qi
            let rhs = gate.mul_add(ctx_gate, r1i_at_gamma, qi_constants[z], rhs);

            // rhs = rhs + r2i(gamma) * cyclo(gamma)
            let rhs = gate.mul_add(ctx_gate, r2i_at_gamma, cyclo_at_gamma_assigned, rhs);
            let lhs = ct0is_at_gamma_assigned[z];

            // LHS(gamma) = RHS(gamma)
            let res = gate.is_equal(ctx_gate, lhs, rhs);
            gate.assert_is_const(ctx_gate, &res, &F::from(1));
        }

        for z in 0..ct0is.len() {
            let pk1_u = gate.mul_add(ctx_gate, pk1_qi_at_gamma[z], u_at_gamma, e1_at_gamma);

            let p1i_at_gamma = p1is_assigned[z].enforce_eval_at_gamma(ctx_rlc, rlc);
            let p2i_at_gamma = p2is_assigned[z].enforce_eval_at_gamma(ctx_rlc, rlc);

            //rhs = pk1_u + p2i * cyclo(gamma)
            let rhs = gate.mul_add(ctx_gate, p2i_at_gamma, cyclo_at_gamma_assigned, pk1_u);

            //rhs = rhs + p1s * qi
            let rhs = gate.mul_add(ctx_gate, p1i_at_gamma, qi_constants[z], rhs);

            let lhs = ct1is_at_gamma_assigned[z];

            let res = gate.is_equal(ctx_gate, lhs, rhs);
            gate.assert_is_const(ctx_gate, &res, &F::from(1));
        }
    }
}
#[cfg(test)]
mod test {

    use std::{fs::File, io::Read};

    use axiom_eth::{
        halo2_proofs::plonk::{keygen_pk, keygen_vk},
        halo2curves::bn256::Fr,
        rlc::{circuit::builder::RlcCircuitBuilder, utils::executor::RlcExecutor},
    };
    use halo2_base::{
        gates::circuit::CircuitBuilderStage,
        utils::{
            fs::gen_srs,
            testing::{check_proof, gen_proof},
        },
    };

    use super::{test_params, BfvPkEncryptionCircuit};

    #[test]
    fn test_pk_enc_valid() {
        let file_path_zeros = "src/data/pk_enc_data/pk_enc_1024_15x60_65537_zeroes.json";
        let mut file = File::open(file_path_zeros).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        let empty_pk_enc_circuit = serde_json::from_str::<BfvPkEncryptionCircuit>(&data).unwrap();

        // 2. Generate (unsafe) trusted setup parameters
        // Here we are setting a small k for optimization purposes
        let k = 14;
        let kzg_params = gen_srs(k as u32);

        // 3. Build the circuit for key generation,
        let mut key_gen_builder =
            RlcCircuitBuilder::<Fr>::from_stage(CircuitBuilderStage::Keygen, 0).use_k(k);
        key_gen_builder.base.set_lookup_bits(k - 1); // lookup bits set to `k-1` as suggested [here](https://docs.axiom.xyz/protocol/zero-knowledge-proofs/getting-started-with-halo2#technical-detail-how-to-choose-lookup_bits)

        let rlc_circuit = RlcExecutor::new(key_gen_builder, empty_pk_enc_circuit.clone());

        // The parameters are auto configured by halo2 lib to fit all the columns into the `k`-sized table
        let rlc_circuit_params = rlc_circuit.0.calculate_params(Some(9));

        // 4. Generate the verification key and the proving key
        let vk = keygen_vk(&kzg_params, &rlc_circuit).unwrap();
        let pk = keygen_pk(&kzg_params, vk, &rlc_circuit).unwrap();
        let break_points = rlc_circuit.0.builder.borrow().break_points();
        drop(rlc_circuit);

        // 5. Generate the proof, here we pass the actual inputs
        let mut proof_gen_builder: RlcCircuitBuilder<Fr> =
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Prover, 0)
                .use_params(rlc_circuit_params);
        proof_gen_builder.base.set_lookup_bits(k - 1);

        let file_path = "src/data/pk_enc_data/pk_enc_1024_15x60_65537.json";
        let mut file = File::open(file_path).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        let sk_enc_circuit = serde_json::from_str::<BfvPkEncryptionCircuit>(&data).unwrap();

        let rlc_circuit = RlcExecutor::new(proof_gen_builder, sk_enc_circuit.clone());

        rlc_circuit
            .0
            .builder
            .borrow_mut()
            .set_break_points(break_points);
        let proof = gen_proof(&kzg_params, &pk, rlc_circuit);

        // 6. Verify the proof
        check_proof(&kzg_params, pk.get_vk(), &proof, true);
    }
}
