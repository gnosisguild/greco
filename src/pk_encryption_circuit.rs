use core::assert;

use axiom_eth::halo2_base::{
    gates::{circuit::BaseCircuitParams, GateInstructions, RangeChip, RangeInstructions},
    utils::ScalarField,
    QuantumCell::Constant,
};
use axiom_eth::rlc::{
    chip::RlcChip,
    circuit::{builder::RlcCircuitBuilder, instructions::RlcCircuitInstructions, RlcCircuitParams},
};
use serde::Deserialize;

use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::str::FromStr;

use crate::{
    constants::pk_enc_constants::pk_enc_constants_2048_1x52_1032193::{
        //constants::pk_enc_constants::pk_enc_constants_1024_15x60_65537::{
        E_BOUND,
        K0IS,
        K1_LOW_BOUND,
        K1_UP_BOUND,
        N,
        P1_BOUNDS,
        P2_BOUNDS,
        PK_BOUND,
        QIS,
        R1_LOW_BOUNDS,
        R1_UP_BOUNDS,
        R2_BOUNDS,
        U_BOUND,
    },
    poly::{Poly, PolyAssigned},
};

/// Helper function to define the parameters of the RlcCircuit. This is a non-optimized configuration that makes use of a single advice column. Use this for testing purposes only.
pub fn test_params() -> RlcCircuitParams {
    RlcCircuitParams {
        base: BaseCircuitParams {
            k: 22,
            num_advice_per_phase: vec![1, 1],
            num_fixed: 1,
            num_lookup_advice_per_phase: vec![0, 1],
            lookup_bits: Some(8),
            num_instance_columns: 1,
        },
        num_rlc_columns: 1,
    }
}
/// `BfvPkEncryptionCircuit` is a circuit that checks the correct formation of a ciphertext resulting from BFV public key encryption
/// All the polynomials coefficients and scalars are normalized to be in the range `[0, p)` where p is the modulus of the prime field of the circuit
///
/// pk_q1 = ( pk0i , pk1i )=( [ai*s + E] , -ai )
/// # Parameters:
/// * `pk0i`: publicly polynomial created by secret polynomial ([ai*s + E] )
/// * `pk1i`: publicly polynomial created by polynomial (-[ai])
/// * `u`: secret polynomial, sampled from ternary distribution.
/// * `e0`: error polynomial, sampled from discrete Gaussian distribution.
/// * `e1`: error polynomial, sampled from discrete Gaussian distribution.
/// * `k1`: scaled message polynomial.
/// * `r2is`: list of r2i polynomials for each i-th CRT basis .
/// * `r1is`: list of r1i polynomials for each CRT i-th CRT basis.
/// * `p2is`: list of p2i polynomials for each i-th CRT basis.
/// * `p1is`: list of p1i polynomials for each i-th CRT basis.
/// * `ct0is`: list of ct0i (first component of the ciphertext cti) polynomials for each CRT i-th CRT basis.
/// * `ct1is`: list of ct1i (second component of the ciphertext cti) polynomials for each CRT i-th CRT basis.

#[derive(Deserialize, Clone)]
pub struct BfvPkEncryptionCircuit {
    pk0i: Vec<Vec<String>>,
    pk1i: Vec<Vec<String>>,
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

impl BfvPkEncryptionCircuit {
    pub fn create_empty_circuit(num_moduli: usize, degree: usize) -> Self {
        let zero_str = String::from("0");

        BfvPkEncryptionCircuit {
            pk0i: vec![vec![zero_str.clone(); degree]; num_moduli],
            pk1i: vec![vec![zero_str.clone(); degree]; num_moduli],
            ct0is: vec![vec![zero_str.clone(); degree]; num_moduli],
            ct1is: vec![vec![zero_str.clone(); degree]; num_moduli],
            r1is: vec![vec![zero_str.clone(); 2 * (degree - 1) + 1]; num_moduli],
            r2is: vec![vec![zero_str.clone(); degree - 1]; num_moduli],
            p1is: vec![vec![zero_str.clone(); 2 * (degree - 1) + 1]; num_moduli],
            p2is: vec![vec![zero_str.clone(); degree - 1]; num_moduli],
            u: vec![zero_str.clone(); degree],
            e0: vec![zero_str.clone(); degree],
            e1: vec![zero_str.clone(); degree],
            k1: vec![zero_str.clone(); degree],
        }
    }
}

/// Payload returned by the first phase of the circuit to be reused in the second phase
pub struct Payload<F: ScalarField> {
    pk0i_assigned: Vec<PolyAssigned<F>>,
    pk1i_assigned: Vec<PolyAssigned<F>>,
    u_assigned: PolyAssigned<F>,
    e0_assigned: PolyAssigned<F>,
    e1_assigned: PolyAssigned<F>,
    k1_assigned: PolyAssigned<F>,
    r2is_assigned: Vec<PolyAssigned<F>>,
    r1is_assigned: Vec<PolyAssigned<F>>,
    p2is_assigned: Vec<PolyAssigned<F>>,
    p1is_assigned: Vec<PolyAssigned<F>>,
    ct0is_assigned: Vec<PolyAssigned<F>>,
    ct1is_assigned: Vec<PolyAssigned<F>>,
}

impl BfvPkEncryptionCircuit {
    /// Helper functions added for checking if the public polynomials (Ct0i, Ct1i, Pk0i, Pk1i) are of the right degree. Also for checking if the coefficients of these polynomials are within the right bounds. Note these checks are necessary for the zkp to be valid. That is, it is true that we can always reduce these polynomials so as to be in R_qi, and hence have the right degree and the coefficients in the right bounds, but accepting polynomials that are not already in R_qi will allow an attacker to forge a valid proof for an invalid ciphertext, exploiting the fact that the equations checked with the zkp live in Z_p and not Z_qi.
    ///
    ///Note as well that the values K0,i have to be the right ones (negation of the inverse of t modulu qi), but we assume that these are calculated by the verifier and not retrieved from the prover, hence no check is performed. For instance, if we have an evm verifier, then this verifier can for instance calculate them himself or retrieve them from the blockchain if these are published.

    pub fn check_polynomial_bounds(&self) {
        //Note we are hardwiring the prime field of the snark in here. There should be a cleaner
        //way to do this.
        let p = BigInt::from_str(
            "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        )
        .expect("Invalid prime number p");

        //Note we are using PK_BOUND for all polynomials we are checking, as all these polynomials
        //have the same bounds.
        let bounds: Vec<BigInt> = PK_BOUND.iter().map(|&b| BigInt::from(b)).collect(); // Convert PK_BOUND to BigInt
        let expected_length = bounds.len();

        // Ensure all polynomials have the same number of rows as PK_BOUND
        if self.pk0i.len() != expected_length
            || self.pk1i.len() != expected_length
            || self.ct0is.len() != expected_length
            || self.ct1is.len() != expected_length
        {
            panic!(
                "Mismatch in polynomial row counts: Expected {}, but got pk0i={}, pk1i={}, ct0is={}, ct1is={}",
                expected_length,
                self.pk0i.len(),
                self.pk1i.len(),
                self.ct0is.len(),
                self.ct1is.len()
            );
        }

        let polynomials = [
            ("pk0i", &self.pk0i),
            ("pk1i", &self.pk1i),
            ("ct0is", &self.ct0is),
            ("ct1is", &self.ct1is),
        ];

        for (name, poly) in polynomials.iter() {
            for (row_idx, row) in poly.iter().enumerate() {
                if row.len() != N {
                    panic!(
                        "Row size mismatch in {}: Expected {}, but row {} has size {}",
                        name,
                        N,
                        row_idx,
                        row.len()
                    );
                }
            }
        }

        fn is_in_bound(value: &str, bound: &BigInt, p: &BigInt) -> bool {
            if let Ok(num) = BigInt::from_str(value) {
                let adjusted = (num + bound) % p; // Apply transformation
                adjusted >= BigInt::zero() && adjusted <= (bound * &BigInt::from(2))
            // Check [0, 2*PK_BOUND]
            } else {
                false // If parsing fails, consider it invalid
            }
        }

        // Helper function to check polynomial bounds per row
        fn check_poly_bounds(
            poly: &Vec<Vec<String>>,
            bounds: &[BigInt],
            p: &BigInt,
        ) -> Vec<(usize, usize, String)> {
            poly.iter()
                .enumerate()
                .flat_map(|(row, row_values)| {
                    if row >= bounds.len() {
                        return vec![]; // Avoid out-of-bounds indexing for PK_BOUND
                    }
                    let bound = &bounds[row];
                    row_values
                        .iter()
                        .enumerate()
                        .filter_map(|(col, value)| {
                            if !is_in_bound(value, bound, p) {
                                Some((row, col, value.clone())) // Collect invalid values
                            } else {
                                None
                            }
                        })
                        .collect::<Vec<_>>()
                })
                .collect()
        }

        // Validate pk0i, pk1i, ct0is, ct1is
        let pk0i_invalid = check_poly_bounds(&self.pk0i, &bounds, &p);
        let pk1i_invalid = check_poly_bounds(&self.pk1i, &bounds, &p);
        let ct0is_invalid = check_poly_bounds(&self.ct0is, &bounds, &p);
        let ct1is_invalid = check_poly_bounds(&self.ct1is, &bounds, &p);

        // Collect all errors
        let mut error_messages = Vec::new();

        if !pk0i_invalid.is_empty() {
            error_messages.push(format!("Invalid pk0i values at {:?}", pk0i_invalid));
        }
        if !pk1i_invalid.is_empty() {
            error_messages.push(format!("Invalid pk1i values at {:?}", pk1i_invalid));
        }
        if !ct0is_invalid.is_empty() {
            error_messages.push(format!("Invalid ct0is values at {:?}", ct0is_invalid));
        }
        if !ct1is_invalid.is_empty() {
            error_messages.push(format!("Invalid ct1is values at {:?}", ct1is_invalid));
        }

        // If there are errors, panic with detailed message
        if !error_messages.is_empty() {
            panic!("{}", error_messages.join(" | "));
        }
    }
}

impl<F: ScalarField> RlcCircuitInstructions<F> for BfvPkEncryptionCircuit {
    type FirstPhasePayload = Payload<F>;

    /// #### Phase 0

    /// In this phase, the polynomials for each matrix $S^j_i$ are assigned to the circuit (j from 0 to 1 refering to two equations, one for ct0 and one for ct1. i from 0 to l-1, where l is the number of qi's). Namely:
    /// * polynomials `u`,'e1, `e0`, `k1` are assigned to the witness table. This has to be done only once as these polynomial are common to each $S_i$ matrix
    /// * polynomials `r1i`, `r2i` are assigned to the witness table for each $S^0_i$ matrix
    /// * polynomial 'ct0is', `pk0i` are assigned to the witness table for each $S^0_i$ matrix and exposed as public inputs
    /// * polynomials `p1i`,`p2i` are assigned to the witness table for each $S^1_i$ matrix
    /// * polynomials 'ct1is`,`pk1i are assigned to the witness table for each $S^1_i$ matrix and exposed as public inputs

    /// Witness values are element of the finite field $\mod{p}$. Negative coefficients $-z$ are assigned as field elements $p - z$.

    /// At the end of phase 0, the witness generated so far is interpolated into a polynomial and committed by the prover. The hash of this commitment is used as challenge and will be used as a source of randomness $\gamma$ in Phase 1. This feature is made available by Halo2 [Challenge API](https://hackmd.io/@axiom/SJw3p-qX3).

    fn virtual_assign_phase0(
        &self,
        builder: &mut RlcCircuitBuilder<F>,
        _: &RangeChip<F>,
    ) -> Self::FirstPhasePayload {
        let ctx = builder.base.main(0);

        let mut public_inputs = vec![];

        let pk0i = self
            .pk0i
            .iter()
            .map(|pk0i| Poly::<F>::new(pk0i.clone()))
            .collect::<Vec<_>>();
        let pk0i_assigned = pk0i
            .into_iter()
            .map(|pk0i| PolyAssigned::new(ctx, pk0i))
            .collect::<Vec<_>>();

        let pk1i = self
            .pk1i
            .iter()
            .map(|pk1i| Poly::<F>::new(pk1i.clone()))
            .collect::<Vec<_>>();
        let pk1i_assigned = pk1i
            .into_iter()
            .map(|pk1i| PolyAssigned::new(ctx, pk1i))
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

        let ct0is_assigned = self
            .ct0is
            .iter()
            .map(|ct0is| {
                let ct0is = Poly::<F>::new(ct0is.clone());
                PolyAssigned::new(ctx, ct0is)
            })
            .collect::<Vec<_>>();

        let ct1is_assigned = self
            .ct1is
            .iter()
            .map(|ct1is| {
                let ct1is = Poly::<F>::new(ct1is.clone());
                PolyAssigned::new(ctx, ct1is)
            })
            .collect::<Vec<_>>();

        for pk0 in pk0i_assigned.iter() {
            for assigned_coefficient in &pk0.assigned_coefficients {
                public_inputs.push(*assigned_coefficient);
            }
        }
        for pk1 in pk1i_assigned.iter() {
            for assigned_coefficient in &pk1.assigned_coefficients {
                public_inputs.push(*assigned_coefficient);
            }
        }
        for ct0 in ct0is_assigned.iter() {
            for assigned_coefficient in &ct0.assigned_coefficients {
                public_inputs.push(*assigned_coefficient);
            }
        }
        for ct1 in ct1is_assigned.iter() {
            for assigned_coefficient in &ct1.assigned_coefficients {
                public_inputs.push(*assigned_coefficient);
            }
        }

        builder.base.assigned_instances[0] = public_inputs;

        Payload {
            pk0i_assigned,
            pk1i_assigned,
            u_assigned,
            e0_assigned,
            e1_assigned,
            k1_assigned,
            r2is_assigned,
            r1is_assigned,
            p2is_assigned,
            p1is_assigned,
            ct0is_assigned,
            ct1is_assigned,
        }
    }

    /// #### Phase 1

    /// In this phase, the following two core constraints are enforced:

    /// - The coefficients of $S^j_i$ are in the expected range.
    /// - $U^j_i(\gamma) \times S^j_i(\gamma) =Ct_{j,i}(\gamma)$

    /// ##### Range Check

    /// The coefficients of the private polynomials from each $i$-th matrix $S^j_i$ are checked to be in the correct range
    /// * Range check polynomials `u`, `e0`,`e1`,`k1`. This has to be done only once as these polynomial are common to each $S^j_i$ matrix
    /// * Range check polynomials `r1i`, `r2i` for each $S^0_i$ matrix
    /// * Range check polynomials `p1i`, `p2i` for each $S^1_i$ matrix

    /// Since negative coefficients `-z` are assigned as `p - z` to the circuit, this might result in very large coefficients. Performing the range check on such large coefficients requires large lookup tables. To avoid this, the coefficients (both negative and positive) are shifted by a constant to make them positive and then perform the range check.

    /// ##### Evaluation at $\gamma$ Constraint

    /// * Constrain the evaluation of the polynomials `u`, `e0`, `e1`, `k1` at $\gamma$. This has to be done only once as these polynomial are common to each $S_i$ matrix
    /// * Constrain the evaluation of the polynomials `r1i`, `r2i` at $\gamma$ for each $S^0_i$ matrix
    /// * Constrain the evaluation of the polynomials `p1i`, `p2i` at $\gamma$ for each $S^1_i$ matrix
    /// Constrain the evaluation of the polynomials `pk0i`, `ct0i` at $\gamma$ for each $U^0_i$ matrix
    /// Constrain the evaluation of the polynomials `pk1i`, `ct1i` at $\gamma$ for each $U^1_i$ matrix
    /// * Constrain the evaluation of the polynomials `cyclo` at $\gamma$ . This has to be done only once as the cyclotomic polynomial is common to each $U_i$ matrix

    /// ##### Correct Encryption Constraint

    /// It is needed to prove that $U^j_i(\gamma) \times S^j_i(\gamma) =Ct_{j,i}(\gamma)$. This can be rewritten as (for the case of j = 0 for instance) `ct0i = ct0i_hat + r1i * qi + r2i * cyclo`, where `ct0i_hat = pk0i * u + e0 + k1 * k0i`.

    /// This constrained is enforced by proving that `LHS(gamma) = RHS(gamma)`. According to the Schwartz-Zippel lemma, if this relation between polynomial when evaluated at a random point holds true, then then the polynomials are identical with high probability. Note that `qi` and `k0i` (for each $U_i$ matrix) are constants to the circuit encoded during key generation.
    /// * Constrain that `ct0i(gamma) = pk0i(gamma) * u(gamma) + e0(gamma) + k1(gamma) * k0i + r1i(gamma) * qi + r2i(gamma) * cyclo(gamma)` for each $i$-th CRT basis
    ///

    fn virtual_assign_phase1(
        builder: &mut RlcCircuitBuilder<F>,
        range: &RangeChip<F>,
        rlc: &RlcChip<F>,
        payload: Self::FirstPhasePayload,
    ) {
        let Payload {
            pk0i_assigned,
            pk1i_assigned,
            u_assigned,
            e0_assigned,
            e1_assigned,
            k1_assigned,
            r2is_assigned,
            r1is_assigned,
            p2is_assigned,
            p1is_assigned,
            ct0is_assigned,
            ct1is_assigned,
        } = payload;

        // ASSIGNMENT

        let (ctx_gate, ctx_rlc) = builder.rlc_ctx_pair();
        let gate = range.gate();

        let mut qi_constants = vec![];
        let mut k0i_constants = vec![];

        for z in 0..ct0is_assigned.len() {
            let qi_constant = Constant(F::from_str_vartime(QIS[z]).unwrap());
            qi_constants.push(qi_constant);

            let k0i_constant = Constant(F::from_str_vartime(K0IS[z]).unwrap());
            k0i_constants.push(k0i_constant);
        }

        // cyclo poly is equal to x^N + 1
        let bits_used = (usize::BITS as usize) - (N.leading_zeros() as usize);
        rlc.load_rlc_cache((ctx_gate, ctx_rlc), gate, bits_used);
        let cyclo_at_gamma_assigned = rlc.rlc_pow_fixed(ctx_gate, gate, N);
        let cyclo_at_gamma_assigned =
            gate.add(ctx_gate, cyclo_at_gamma_assigned, Constant(F::from(1)));

        u_assigned.range_check_1bound(ctx_gate, range, U_BOUND);
        e0_assigned.range_check_1bound(ctx_gate, range, E_BOUND);
        e1_assigned.range_check_1bound(ctx_gate, range, E_BOUND);
        k1_assigned.range_check_2bounds(ctx_gate, range, K1_LOW_BOUND, K1_UP_BOUND);

        let _ = pk0i_assigned
            .iter()
            .enumerate()
            .map(|(i, pk_assigned)| pk_assigned.range_check_1bound(ctx_gate, range, PK_BOUND[i]));

        let _ = pk1i_assigned
            .iter()
            .enumerate()
            .map(|(i, pk_assigned)| pk_assigned.range_check_1bound(ctx_gate, range, PK_BOUND[i]));

        for z in 0..ct0is_assigned.len() {
            r2is_assigned[z].range_check_1bound(ctx_gate, range, R2_BOUNDS[z]);
            r1is_assigned[z].range_check_2bounds(
                ctx_gate,
                range,
                R1_LOW_BOUNDS[z],
                R1_UP_BOUNDS[z],
            );
            p2is_assigned[z].range_check_1bound(ctx_gate, range, P2_BOUNDS[z]);
            p1is_assigned[z].range_check_1bound(ctx_gate, range, P1_BOUNDS[z]);
        }

        let u_at_gamma = u_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);
        let e0_at_gamma = e0_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);
        let e1_at_gamma = e1_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);
        let k1_at_gamma = k1_assigned.enforce_eval_at_gamma(ctx_rlc, rlc);
        let pk0i_at_gamma = pk0i_assigned
            .iter()
            .map(|pk_assigned| pk_assigned.enforce_eval_at_gamma(ctx_rlc, rlc))
            .collect::<Vec<_>>();
        let pk1i_at_gamma = pk1i_assigned
            .iter()
            .map(|pk_assigned| pk_assigned.enforce_eval_at_gamma(ctx_rlc, rlc))
            .collect::<Vec<_>>();
        let gate = range.gate();

        // For each `i` Prove that LHS(gamma) = RHS(gamma)
        // pk0_u = pk0i(gamma) * u(gamma) + e0(gamma)
        // LHS = ct0i(gamma)
        // RHS = pk0_u  + k1(gamma) * k0i + r1i(gamma) * qi + r2i(gamma) * cyclo(gamma)

        for z in 0..ct0is_assigned.len() {
            let pk0_u = gate.mul_add(ctx_gate, pk0i_at_gamma[z], u_at_gamma, e0_at_gamma);
            let r1i_at_gamma = r1is_assigned[z].enforce_eval_at_gamma(ctx_rlc, rlc);
            let r2i_at_gamma = r2is_assigned[z].enforce_eval_at_gamma(ctx_rlc, rlc);

            // CORRECT ENCRYPTION CONSTRAINT

            // rhs = pk0_u + k1(gamma) * k0i
            let rhs = gate.mul_add(ctx_gate, k1_at_gamma, k0i_constants[z], pk0_u);

            // rhs = rhs + r1i(gamma) * qi
            let rhs = gate.mul_add(ctx_gate, r1i_at_gamma, qi_constants[z], rhs);

            // rhs = rhs + r2i(gamma) * cyclo(gamma)
            let rhs = gate.mul_add(ctx_gate, r2i_at_gamma, cyclo_at_gamma_assigned, rhs);
            let lhs = ct0is_assigned[z].enforce_eval_at_gamma(ctx_rlc, rlc);

            // LHS(gamma) = RHS(gamma)
            let res = gate.is_equal(ctx_gate, lhs, rhs);
            gate.assert_is_const(ctx_gate, &res, &F::from(1));
        }

        for z in 0..ct1is_assigned.len() {
            let pk1_u = gate.mul_add(ctx_gate, pk1i_at_gamma[z], u_at_gamma, e1_at_gamma);

            let p1i_at_gamma = p1is_assigned[z].enforce_eval_at_gamma(ctx_rlc, rlc);
            let p2i_at_gamma = p2is_assigned[z].enforce_eval_at_gamma(ctx_rlc, rlc);

            //rhs = pk1_u + p2i * cyclo(gamma)
            let rhs = gate.mul_add(ctx_gate, p2i_at_gamma, cyclo_at_gamma_assigned, pk1_u);

            //rhs = rhs + p1s * qi
            let rhs = gate.mul_add(ctx_gate, p1i_at_gamma, qi_constants[z], rhs);

            let lhs = ct1is_assigned[z].enforce_eval_at_gamma(ctx_rlc, rlc);

            let res = gate.is_equal(ctx_gate, lhs, rhs);
            gate.assert_is_const(ctx_gate, &res, &F::from(1));
        }
    }

    fn instances(&self) -> Vec<Vec<F>> {
        let mut instance = vec![];
        for pk0 in self.pk0i.iter() {
            let pk0_poly = Poly::<F>::new(pk0.clone());
            instance.extend(pk0_poly.coefficients);
        }
        for pk1 in self.pk1i.iter() {
            let pk1_poly = Poly::<F>::new(pk1.clone());
            instance.extend(pk1_poly.coefficients);
        }
        for ct0i in self.ct0is.iter() {
            let ct0i_poly = Poly::<F>::new(ct0i.clone());
            instance.extend(ct0i_poly.coefficients);
        }
        for ct1i in self.ct1is.iter() {
            let ct1i_poly = Poly::<F>::new(ct1i.clone());
            instance.extend(ct1i_poly.coefficients);
        }
        vec![instance]
    }
}
#[cfg(test)]
mod test {
    use std::fs::File;
    use std::io::Read;
    use std::io::Write;

    use axiom_eth::halo2_base::{
        gates::circuit::CircuitBuilderStage,
        halo2_proofs::{
            dev::{FailureLocation, MockProver, VerifyFailure},
            halo2curves::bn256::{Bn256, Fr},
            plonk::{create_proof, keygen_pk, keygen_vk, verify_proof, Any, SecondPhase},
            poly::kzg::{
                commitment::ParamsKZG,
                multiopen::{ProverSHPLONK, VerifierSHPLONK},
                strategy::SingleStrategy,
            },
            transcript::TranscriptWriterBuffer,
        },
        utils::{
            fs::gen_srs,
            testing::{check_proof_with_instances, gen_proof_with_instances},
        },
    };
    use halo2_solidity_verifier::{
        compile_solidity, encode_calldata, BatchOpenScheme::Bdfg21, Evm, Keccak256Transcript,
        SolidityGenerator,
    };
    use rand::{rngs::OsRng, rngs::StdRng, RngCore, SeedableRng};

    use axiom_eth::rlc::{
        circuit::{builder::RlcCircuitBuilder, instructions::RlcCircuitInstructions},
        utils::executor::RlcExecutor,
    };

    //use crate::constants::pk_enc_constants::pk_enc_constants_1024_15x60_65537::R1_LOW_BOUNDS;

    use super::{test_params, BfvPkEncryptionCircuit};

    #[test]
    fn test_pk_enc_valid() {
        let file_path_zeroes = "src/data/pk_enc_data/pk_enc_1024_2x52_2048_zeroes.json";
        //let file_path_zeroes = "src/data/pk_enc_data/pk_enc_1024_15x60_65537_zeroes.json";
        let mut file = File::open(file_path_zeroes).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        let empty_pk_enc_circuit = serde_json::from_str::<BfvPkEncryptionCircuit>(&data).unwrap();

        // 2. Generate (unsafe) trusted setup parameters
        // Here we are setting a small k for optimization purposes
        let k = 16;
        let kzg_params = gen_srs(k as u32);

        // 3. Build the circuit for key generation,
        let mut key_gen_builder =
            RlcCircuitBuilder::<Fr>::from_stage(CircuitBuilderStage::Keygen, 0).use_k(k);
        key_gen_builder.base.set_lookup_bits(k - 1); // lookup bits set to `k-1` as suggested [here](https://docs.axiom.xyz/protocol/zero-knowledge-proofs/getting-started-with-halo2#technical-detail-how-to-choose-lookup_bits)
        key_gen_builder.base.set_instance_columns(1);
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
            RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0)
                .use_params(rlc_circuit_params);
        proof_gen_builder.base.set_lookup_bits(k - 1);
        proof_gen_builder.base.set_instance_columns(1);

        let file_path = "src/data/pk_enc_data/pk_enc_1024_2x52_2048.json";
        //let file_path = "src/data/pk_enc_data/pk_enc_1024_15x60_65537.json";
        let mut file = File::open(file_path).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();
        let pk_enc_circuit = serde_json::from_str::<BfvPkEncryptionCircuit>(&data).unwrap();

        let rlc_circuit = RlcExecutor::new(proof_gen_builder, pk_enc_circuit.clone());

        rlc_circuit
            .0
            .builder
            .borrow_mut()
            .set_break_points(break_points);
        let instances = pk_enc_circuit.instances();
        let proof = gen_proof_with_instances(&kzg_params, &pk, rlc_circuit, &[&instances[0]]);

        // 6. Verify the proof. Note the first check is as well part of verifying the zkp

        pk_enc_circuit.check_polynomial_bounds();
        check_proof_with_instances(&kzg_params, pk.get_vk(), &proof, &[&instances[0]], true);
    }

    #[test]
    pub fn test_pk_enc_full_prover() {
        // --------------------------------------------------
        // (A) Generate a proof & verify it locally
        // --------------------------------------------------
        // Zero file for keygen circuit sizing
        let empty_pk_enc_circuit = BfvPkEncryptionCircuit::create_empty_circuit(1, 2048);

        let k = 15;
        let kzg_params = gen_srs(k);

        // Build an RLC circuit for KeyGen
        let mut key_gen_builder =
            RlcCircuitBuilder::<Fr>::from_stage(CircuitBuilderStage::Keygen, 0).use_k(k as usize);
        key_gen_builder.base.set_lookup_bits((k - 1) as usize);
        key_gen_builder.base.set_instance_columns(1);

        let rlc_circuit_for_keygen =
            RlcExecutor::new(key_gen_builder, empty_pk_enc_circuit.clone());
        let rlc_circuit_params = rlc_circuit_for_keygen.0.calculate_params(Some(9));

        // Keygen VerifyingKey / ProvingKey
        let vk = keygen_vk(&kzg_params, &rlc_circuit_for_keygen).unwrap();
        let pk = keygen_pk(&kzg_params, vk, &rlc_circuit_for_keygen).unwrap();
        let actual_num_instance_columns = pk.get_vk().cs().num_instance_columns();
        println!("VerifyingKey says num_instance_columns = {actual_num_instance_columns}");

        let break_points = rlc_circuit_for_keygen.0.builder.borrow().break_points();
        drop(rlc_circuit_for_keygen);

        // Load the real data from JSON
        let file_path = "src/data/pk_enc_data/pk_enc_2048_1x52_1032193.json";
        let pk_enc_circuit: BfvPkEncryptionCircuit =
            serde_json::from_reader(File::open(file_path).unwrap()).unwrap();
        let instances: Vec<Vec<Fr>> = pk_enc_circuit.instances();

        println!("instances.len() = {}", instances.len());
        println!("instances[0].len() = {}", instances[0].len());

        // Build the RLC circuit for the real data
        let mut builder = RlcCircuitBuilder::from_stage(CircuitBuilderStage::Prover, 0)
            .use_params(rlc_circuit_params.clone());
        // lookup bits set to `k-1` as suggested [here](https://docs.axiom.xyz/protocol/zero-knowledge-proofs/getting-started-with-halo2#technical-detail-how-to-choose-lookup_bits)
        builder.base.set_lookup_bits((k - 1) as usize);
        builder.base.set_instance_columns(1);

        let rlc_prover_circuit = RlcExecutor::new(builder, pk_enc_circuit.clone());
        rlc_prover_circuit
            .0
            .builder
            .borrow_mut()
            .set_break_points(break_points);

        // Create a proof
        let mut rng = StdRng::seed_from_u64(OsRng.next_u64());
        let instance_refs = vec![instances[0].as_slice()];

        let proof = {
            let mut transcript = Keccak256Transcript::new(Vec::new());
            create_proof::<_, ProverSHPLONK<_>, _, _, _, _>(
                &kzg_params,
                &pk,
                &[rlc_prover_circuit],
                &[&instance_refs],
                &mut rng,
                &mut transcript,
            )
            .unwrap();
            transcript.finalize()
        };

        println!("E2E proof size = {} bytes", proof.len());

        // Verify it locally
        let result = {
            let mut transcript = Keccak256Transcript::new(proof.as_slice());
            verify_proof::<_, VerifierSHPLONK<_>, _, _, SingleStrategy<_>>(
                &kzg_params,
                pk.get_vk(),
                SingleStrategy::new(&kzg_params),
                &[&instance_refs],
                &mut transcript,
            )
        };
        assert!(result.is_ok());
        println!("Local verification succeeded!");

        // --------------------------------------------------
        // (B) Generate a Solidity verifier & test in an EVM
        // --------------------------------------------------
        let num_public_inputs = instances[0].len();

        let generator = SolidityGenerator::new(&kzg_params, pk.get_vk(), Bdfg21, num_public_inputs);

        // Render the Solidity code as a single contract
        let verifier_solidity: String = generator.render().expect("render contract");

        // Write it to a file
        // let mut file = File::create("./contracts/BfvPKEncryptionVerifier.sol").unwrap();
        // file.write_all(verifier_solidity.as_bytes()).unwrap();
        // println!("Solidity verifier contract written to PKVerifier.sol");

        // Compile it
        let creation_code = compile_solidity(&verifier_solidity);
        let code_size = creation_code.len();
        println!("Verifier creation code size: {}", code_size);

        // Deploy it to a local EVM
        let mut evm = Evm::default();
        let verifier_address = evm.create(creation_code);
        println!("verifier_address = {:?}", verifier_address);

        // Encode the calldata: we have None for "inlined" verifying key in the same contract
        let calldata = encode_calldata(None, &proof, &instances[0]);

        // Call the contract
        let (gas_cost, output) = evm.call(verifier_address, calldata);

        assert_eq!(output.last(), Some(&1u8), "EVM returned 'false'");
        println!("EVM verification success with gas cost = {gas_cost}");
    }
    // #[test]
    // pub fn test_pk_enc_invalid_range() {
    //     // 1. Define the inputs of the circuit

    //     let file_path = "src/data/pk_enc_data/pk_enc_1024_15x60_65537.json";
    //     let mut file = File::open(file_path).unwrap();
    //     let mut data = String::new();
    //     file.read_to_string(&mut data).unwrap();
    //     let mut pk_enc_circuit = serde_json::from_str::<BfvPkEncryptionCircuit>(&data).unwrap();
    //     let instances = pk_enc_circuit.instances();

    //     // 2. Invalidate the circuit by setting the value of a coefficient of the polynomial `r1is[0]` to be out of range
    //     let out_of_range_coeff = R1_LOW_BOUNDS[0] + 1;
    //     pk_enc_circuit.r1is[0][0] = out_of_range_coeff.to_string();

    //     // 3. Build the circuit for MockProver
    //     let rlc_circuit_params = test_params();
    //     let mut mock_builder: RlcCircuitBuilder<Fr> =
    //         RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0)
    //             .use_params(rlc_circuit_params.clone());
    //     mock_builder.base.set_lookup_bits(8); // Set the lookup bits to 8
    //     mock_builder.base.set_instance_columns(1);

    //     let rlc_circuit = RlcExecutor::new(mock_builder, pk_enc_circuit);

    //     // 4. Run the mock prover
    //     let invalid_mock_prover = MockProver::run(
    //         rlc_circuit_params.base.k.try_into().unwrap(),
    //         &rlc_circuit,
    //         instances,
    //     )
    //     .unwrap();

    //     // 5. Assert that the circuit is not satisfied
    //     // In particular, it should fail the range check enforced in the second phase for the first coefficient of r1is[0] and the equality check in the second phase for the 0-th basis
    //     assert_eq!(
    //         invalid_mock_prover.verify(),
    //         Err(vec![
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 115709 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 115719 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914202 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914222 }
    //             },
    //         ])
    //     );
    // }

    // #[test]
    // pub fn test_pk_enc_invalid_polys() {
    //     // 1. Define the inputs of the circuit

    //     let file_path = "src/data/pk_enc_data/pk_enc_1024_15x60_65537.json";
    //     let mut file = File::open(file_path).unwrap();
    //     let mut data = String::new();
    //     file.read_to_string(&mut data).unwrap();
    //     let mut pk_enc_circuit = serde_json::from_str::<BfvPkEncryptionCircuit>(&data).unwrap();

    //     // 2. Invalidate the circuit by setting a different `s` polynomial

    //     let invalid_u = vec!["1".to_string(); 1024];
    //     pk_enc_circuit.u = invalid_u;

    //     // 3. Build the circuit for MockProver
    //     let rlc_circuit_params = test_params();
    //     let mut mock_builder: RlcCircuitBuilder<Fr> =
    //         RlcCircuitBuilder::from_stage(CircuitBuilderStage::Mock, 0)
    //             .use_params(rlc_circuit_params.clone());
    //     mock_builder.base.set_lookup_bits(8); // Set the lookup bits to 8

    //     let instances = pk_enc_circuit.instances();
    //     let rlc_circuit = RlcExecutor::new(mock_builder, pk_enc_circuit);

    //     // 4. Run the mock prover
    //     let invalid_mock_prover = MockProver::run(
    //         rlc_circuit_params.base.k.try_into().unwrap(),
    //         &rlc_circuit,
    //         instances,
    //     )
    //     .unwrap();

    //     // 5. Assert that the circuit is not satisfied
    //     // In particular, it should fail the equality check (LHS=RHS) in the second phase for each i-th CRT basis
    //     assert_eq!(
    //         invalid_mock_prover.verify(),
    //         Err(vec![
    //             VerifyFailure::Permutation {
    //                 column: (Any::Fixed, 1).into(),
    //                 location: FailureLocation::InRegion {
    //                     region: (2, "base+rlc phase 1").into(),
    //                     offset: 1
    //                 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914202 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914222 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914230 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914250 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914258 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914278 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914286 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914306 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914314 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914334 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914342 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914362 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914370 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914390 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914398 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914418 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914426 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914446 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914454 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914474 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914482 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914502 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914510 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914530 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914538 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914558 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914566 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914586 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914594 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914610 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914618 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914634 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914642 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914658 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914666 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914682 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914690 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914706 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914714 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914730 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914738 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914754 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914762 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914778 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914786 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914802 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914810 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914826 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914834 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914850 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914858 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914874 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914882 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914898 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914906 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914922 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914930 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914946 }
    //             },
    //             VerifyFailure::Permutation {
    //                 column: (Any::advice_in(SecondPhase), 1).into(),
    //                 location: FailureLocation::OutsideRegion { row: 2914954 }
    //             }
    //         ])
    //     );
    // }

    #[test]
    #[cfg(feature = "bench")]
    pub fn bench_pk_enc_full_prover() {
        //let file_path = "src/data/pk_enc_data/pk_enc_1024_15x60_65537.json";
        let file_path = "src/data/pk_enc_data/pk_enc_1024_2x52_2048.json";

        pub struct Config {
            kzg_params: ParamsKZG<Bn256>,
            k: usize,
        }

        // Generate unsafe parameters for different values of k
        let mut configs = vec![];
        for k in 16..=21 {
            let kzg_params = gen_srs(k as u32);
            let config = Config { kzg_params, k };
            configs.push(config);
        }

        // Prepare a table to display results
        let mut table = Table::new();
        table.add_row(row![
            "K",
            "VK Generation Time",
            "PK Generation Time",
            "Proof Generation Time",
            "Proof Verification Time"
        ]);

        for config in &configs {
            println!("Running bench for k={}", config.k);
            // 1. Define the inputs of the circuit.
            // Since we are going to use this circuit instance for key gen, we can use an input file in which all the coefficients are set to 0
            //let file_path_zeroes = "src/data/pk_enc_data/pk_enc_1024_15x60_65537_zeroes.json";
            let file_path_zeroes = "src/data/pk_enc_data/pk_enc_1024_2x52_2048_zeroes.json";
            let mut file = File::open(file_path_zeroes).unwrap();
            let mut data = String::new();
            file.read_to_string(&mut data).unwrap();
            let empty_pk_enc_circuit =
                serde_json::from_str::<BfvPkEncryptionCircuit>(&data).unwrap();

            // 2. Build the circuit for key generation,
            let mut key_gen_builder =
                RlcCircuitBuilder::from_stage(CircuitBuilderStage::Keygen, 0).use_k(config.k);
            key_gen_builder.base.set_lookup_bits(config.k - 1); // lookup bits set to `k-1` as suggested [here](https://docs.axiom.xyz/protocol/zero-knowledge-proofs/getting-started-with-halo2#technical-detail-how-to-choose-lookup_bits)

            key_gen_builder.base.set_instance_columns(1);

            let rlc_circuit = RlcExecutor::new(key_gen_builder, empty_pk_enc_circuit.clone());

            // The parameters are auto configured by halo2 lib to fit all the columns into the `k`-sized table
            let rlc_circuit_params = rlc_circuit.0.calculate_params(Some(9));

            // 3. Generate the verification key and the proving key
            let timer = std::time::Instant::now();
            let vk = keygen_vk(&config.kzg_params, &rlc_circuit).unwrap();
            let vk_gen_time = timer.elapsed();
            let timer = std::time::Instant::now();
            let pk = keygen_pk(&config.kzg_params, vk, &rlc_circuit).unwrap();
            let pk_gen_time = timer.elapsed();
            let break_points = rlc_circuit.0.builder.borrow().break_points();
            drop(rlc_circuit);

            // 4. Generate the proof, here we pass the actual inputs
            let mut proof_gen_builder: RlcCircuitBuilder<Fr> =
                RlcCircuitBuilder::from_stage(CircuitBuilderStage::Prover, 0)
                    .use_params(rlc_circuit_params);
            proof_gen_builder.base.set_lookup_bits(config.k - 1);

            //let file_path = "src/data/pk_enc_data/pk_enc_1024_15x60_65537.json";
            let file_path = "src/data/pk_enc_data/pk_enc_1024_2x52_2048.json";
            let mut file = File::open(file_path).unwrap();
            let mut data = String::new();
            file.read_to_string(&mut data).unwrap();
            let pk_enc_circuit = serde_json::from_str::<BfvPkEncryptionCircuit>(&data).unwrap();

            let rlc_circuit = RlcExecutor::new(proof_gen_builder, pk_enc_circuit.clone());

            rlc_circuit
                .0
                .builder
                .borrow_mut()
                .set_break_points(break_points);

            let instances = pk_enc_circuit.instances();

            let timer = std::time::Instant::now();
            let proof =
                gen_proof_with_instances(&config.kzg_params, &pk, rlc_circuit, &[&instances[0]]);
            let proof_gen_time = timer.elapsed();

            // 6. Verify the proof. Note the first check is as well part of verifying the zkp
            let timer = std::time::Instant::now();
            pk_enc_circuit.check_polynomial_bounds();
            check_proof_with_instances(
                &config.kzg_params,
                pk.get_vk(),
                &proof,
                &[&instances[0]],
                true,
            );
            let proof_verification_time = timer.elapsed();

            table.add_row(row![
                config.k,
                format!("{:?}", vk_gen_time),
                format!("{:?}", pk_gen_time),
                format!("{:?}", proof_gen_time),
                format!("{:?}", proof_verification_time)
            ]);
        }
        println!("bfv params: {:?}", file_path);
        table.printstd();
    }
}
