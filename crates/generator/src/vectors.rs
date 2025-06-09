//! Input validation vectors for Greco zero-knowledge proofs.
//!
//! This module contains the core data structure and computation logic for generating
//! input validation vectors required for proving correct BFV encryption in zero-knowledge.

use fhe::bfv::{BfvParameters, Ciphertext, Plaintext, PublicKey};
use fhe_math::{
    rq::{Poly, Representation},
    zq::Modulus,
};
use fhe_traits::*;
use itertools::izip;
use num_bigint::BigInt;
use num_traits::{ToPrimitive, Zero};
use polynomial::*;
use rayon::iter::{ParallelBridge, ParallelIterator};
use serde_json::json;
use std::ops::Deref;
use std::sync::Arc;

use crate::utils::{to_string_1d_vec, to_string_2d_vec};

/// Set of vectors for input validation of a ciphertext
#[derive(Clone, Debug)]
pub struct InputValidationVectors {
    pub pk0is: Vec<Vec<BigInt>>,
    pub pk1is: Vec<Vec<BigInt>>,
    pub ct0is: Vec<Vec<BigInt>>,
    pub ct1is: Vec<Vec<BigInt>>,
    pub r1is: Vec<Vec<BigInt>>,
    pub r2is: Vec<Vec<BigInt>>,
    pub p1is: Vec<Vec<BigInt>>,
    pub p2is: Vec<Vec<BigInt>>,
    pub k0is: Vec<BigInt>,
    pub u: Vec<BigInt>,
    pub e0: Vec<BigInt>,
    pub e1: Vec<BigInt>,
    pub k1: Vec<BigInt>,
}

impl InputValidationVectors {
    /// Create a new `InputValidationVectors` with the given number of moduli and degree.
    ///
    /// # Arguments
    ///
    /// * `num_moduli` - The number of moduli, which determines the number of inner vectors in 2D vectors.
    /// * `degree` - The size of each inner vector in the 2D vectors.
    ///
    /// # Returns
    ///
    /// Returns a new instance of `InputValidationVectors` with all fields initialized to zero.
    pub fn new(num_moduli: usize, degree: usize) -> Self {
        InputValidationVectors {
            pk0is: vec![vec![BigInt::zero(); degree]; num_moduli],
            pk1is: vec![vec![BigInt::zero(); degree]; num_moduli],
            ct0is: vec![vec![BigInt::zero(); degree]; num_moduli],
            ct1is: vec![vec![BigInt::zero(); degree]; num_moduli],
            r1is: vec![vec![BigInt::zero(); 2 * (degree - 1) + 1]; num_moduli],
            r2is: vec![vec![BigInt::zero(); degree - 1]; num_moduli],
            p1is: vec![vec![BigInt::zero(); 2 * (degree - 1) + 1]; num_moduli],
            p2is: vec![vec![BigInt::zero(); degree - 1]; num_moduli],
            k0is: vec![BigInt::zero(); num_moduli],
            u: vec![BigInt::zero(); degree],
            e0: vec![BigInt::zero(); degree],
            e1: vec![BigInt::zero(); degree],
            k1: vec![BigInt::zero(); degree],
        }
    }

    /// Assign and return all of the centered input validation vectors to the ZKP modulus `p`.
    ///
    /// # Arguments
    ///
    /// * `p` - ZKP modulus
    ///
    /// # Returns
    ///
    /// Returns a new `InputValidationVectors` struct with all coefficients reduced modulo `p`.
    pub fn standard_form(&self, p: &BigInt) -> Self {
        InputValidationVectors {
            pk0is: reduce_coefficients_2d(&self.pk0is, p),
            pk1is: reduce_coefficients_2d(&self.pk1is, p),
            ct0is: reduce_coefficients_2d(&self.ct0is, p),
            ct1is: reduce_coefficients_2d(&self.ct1is, p),
            r1is: reduce_coefficients_2d(&self.r1is, p),
            r2is: reduce_coefficients_2d(&self.r2is, p),
            p1is: reduce_coefficients_2d(&self.p1is, p),
            p2is: reduce_coefficients_2d(&self.p2is, p),
            k0is: self.k0is.clone(),
            u: reduce_coefficients(&self.u, p),
            e0: reduce_coefficients(&self.e0, p),
            e1: reduce_coefficients(&self.e1, p),
            k1: reduce_coefficients(&self.k1, p),
        }
    }

    /// Convert the `InputValidationVectors` to a JSON object.
    ///
    /// # Returns
    ///
    /// Returns a `serde_json::Value` representing the JSON serialization of the `InputValidationVectors`.
    pub fn to_json(&self) -> serde_json::Value {
        json!({
            "pk0i": to_string_2d_vec(&self.pk0is),
            "pk1i": to_string_2d_vec(&self.pk1is),
            "u": to_string_1d_vec(&self.u),
            "e0": to_string_1d_vec(&self.e0),
            "e1": to_string_1d_vec(&self.e1),
            "k1": to_string_1d_vec(&self.k1),
            "r2is": to_string_2d_vec(&self.r2is),
            "r1is": to_string_2d_vec(&self.r1is),
            "p2is": to_string_2d_vec(&self.p2is),
            "p1is": to_string_2d_vec(&self.p1is),
            "ct0is": to_string_2d_vec(&self.ct0is),
            "ct1is": to_string_2d_vec(&self.ct1is),
        })
    }

    /// Check whether all members of `self` have the correct length based on the provided `degree` and `num_moduli`.
    ///
    /// # Arguments
    ///
    /// * `num_moduli` - The expected number of moduli (outer vector length).
    /// * `degree` - The expected degree (inner vector length).
    ///
    /// # Returns
    ///
    /// Returns `true` if all vectors have the correct lengths, `false` otherwise.
    pub fn check_correct_lengths(&self, num_moduli: usize, degree: usize) -> bool {
        // Helper function to check 2D vector lengths
        let check_2d_lengths =
            |vec: &Vec<Vec<BigInt>>, expected_outer_len: usize, expected_inner_len: usize| {
                vec.len() == expected_outer_len && vec.iter().all(|v| v.len() == expected_inner_len)
            };

        // Helper function to check 1D vector lengths
        let check_1d_lengths = |vec: &Vec<BigInt>, expected_len: usize| vec.len() == expected_len;

        // Use all to combine all checks into a single statement
        [
            // 2D vector checks
            check_2d_lengths(&self.pk0is, num_moduli, degree),
            check_2d_lengths(&self.pk1is, num_moduli, degree),
            check_2d_lengths(&self.ct0is, num_moduli, degree),
            check_2d_lengths(&self.ct1is, num_moduli, degree),
            check_2d_lengths(&self.r1is, num_moduli, 2 * (degree - 1) + 1),
            check_2d_lengths(&self.r2is, num_moduli, degree - 1),
            check_2d_lengths(&self.p1is, num_moduli, 2 * (degree - 1) + 1),
            check_2d_lengths(&self.p2is, num_moduli, degree - 1),
            // 1D vector checks
            check_1d_lengths(&self.k0is, num_moduli),
            check_1d_lengths(&self.u, degree),
            check_1d_lengths(&self.e0, degree),
            check_1d_lengths(&self.e1, degree),
            check_1d_lengths(&self.k1, degree),
        ]
        .iter()
        .all(|&check| check)
    }

    /// Create the centered validation vectors necessary for creating an input validation proof according to Greco.
    /// For more information, please see https://eprint.iacr.org/2024/594.
    ///
    /// # Arguments
    ///
    /// * `pt` - Plaintext from fhe.rs.
    /// * `u_rns` - Private polynomial used in ciphertext sampled from secret key distribution.
    /// * `e0_rns` - Error polynomial used in ciphertext sampled from error distribution.
    /// * `e1_rns` - Error polynomioal used in cihpertext sampled from error distribution.
    /// * `ct` - Ciphertext from fhe.rs.
    /// * `pk` - Public Key from fhe.rs.
    pub fn compute(
        pt: &Plaintext,
        u_rns: &Poly,
        e0_rns: &Poly,
        e1_rns: &Poly,
        ct: &Ciphertext,
        pk: &PublicKey,
        params: &Arc<BfvParameters>,
    ) -> Result<InputValidationVectors, Box<dyn std::error::Error>> {
        // Get context, plaintext modulus, and degree
        let ctx = params.ctx_at_level(pt.level())?;
        let t = Modulus::new(params.plaintext())?;
        let n: u64 = ctx.degree as u64;

        // Calculate k1 (independent of qi), center and reverse
        let q_mod_t = (ctx.modulus() % t.modulus())
            .to_u64()
            .ok_or_else(|| "Cannot convert BigInt to u64.".to_string())?; // [q]_t
        let mut k1_u64 = pt.value.deref().to_vec(); // m
        t.scalar_mul_vec(&mut k1_u64, q_mod_t); // k1 = [q*m]_t

        let mut k1: Vec<BigInt> = k1_u64.iter().map(|&x| BigInt::from(x)).rev().collect();
        reduce_and_center_coefficients_mut(&mut k1, &BigInt::from(t.modulus()));

        // Extract single vectors of u, e1, and e2 as Vec<BigInt>, center and reverse
        let mut u_rns_copy = u_rns.clone();
        let mut e0_rns_copy = e0_rns.clone();
        let mut e1_rns_copy = e1_rns.clone();

        u_rns_copy.change_representation(Representation::PowerBasis);
        e0_rns_copy.change_representation(Representation::PowerBasis);
        e1_rns_copy.change_representation(Representation::PowerBasis);

        let u: Vec<BigInt> = unsafe {
            ctx.moduli_operators()[0]
                .center_vec_vt(
                    u_rns_copy
                        .coefficients()
                        .row(0)
                        .as_slice()
                        .ok_or_else(|| "Cannot center coefficients.".to_string())?,
                )
                .iter()
                .rev()
                .map(|&x| BigInt::from(x))
                .collect()
        };

        let e0: Vec<BigInt> = unsafe {
            ctx.moduli_operators()[0]
                .center_vec_vt(
                    e0_rns_copy
                        .coefficients()
                        .row(0)
                        .as_slice()
                        .ok_or_else(|| "Cannot center coefficients.".to_string())?,
                )
                .iter()
                .rev()
                .map(|&x| BigInt::from(x))
                .collect()
        };

        let e1: Vec<BigInt> = unsafe {
            ctx.moduli_operators()[0]
                .center_vec_vt(
                    e1_rns_copy
                        .coefficients()
                        .row(0)
                        .as_slice()
                        .ok_or_else(|| "Cannot center coefficients.".to_string())?,
                )
                .iter()
                .rev()
                .map(|&x| BigInt::from(x))
                .collect()
        };

        // Extract and convert ciphertext and plaintext polynomials
        let mut ct0 = ct.c[0].clone();
        let mut ct1 = ct.c[1].clone();
        ct0.change_representation(Representation::PowerBasis);
        ct1.change_representation(Representation::PowerBasis);

        let mut pk0: Poly = pk.c.c[0].clone();
        let mut pk1: Poly = pk.c.c[1].clone();
        pk0.change_representation(Representation::PowerBasis);
        pk1.change_representation(Representation::PowerBasis);

        // Create cyclotomic polynomial x^N + 1
        let mut cyclo = vec![BigInt::from(0u64); (n + 1) as usize];

        cyclo[0] = BigInt::from(1u64); // x^N term
        cyclo[n as usize] = BigInt::from(1u64); // x^0 term

        // Initialize matrices to store results
        let num_moduli = ctx.moduli().len();
        let mut res = InputValidationVectors::new(num_moduli, n as usize);

        let ct0_coeffs = ct0.coefficients();
        let ct1_coeffs = ct1.coefficients();
        let pk0_coeffs = pk0.coefficients();
        let pk1_coeffs = pk1.coefficients();

        let ct0_coeffs_rows = ct0_coeffs.rows();
        let ct1_coeffs_rows = ct1_coeffs.rows();
        let pk0_coeffs_rows = pk0_coeffs.rows();
        let pk1_coeffs_rows = pk1_coeffs.rows();

        // Perform the main computation logic
        let results: Vec<(
            usize,
            Vec<BigInt>,
            Vec<BigInt>,
            BigInt,
            Vec<BigInt>,
            Vec<BigInt>,
            Vec<BigInt>,
            Vec<BigInt>,
            Vec<BigInt>,
            Vec<BigInt>,
        )> = izip!(
            ctx.moduli_operators(),
            ct0_coeffs_rows,
            ct1_coeffs_rows,
            pk0_coeffs_rows,
            pk1_coeffs_rows,
        )
        .enumerate()
        .par_bridge()
        .map(
            |(i, (qi, ct0_coeffs, ct1_coeffs, pk0_coeffs, pk1_coeffs))| {
                // --------------------------------------------------- ct0i ---------------------------------------------------

                // Convert to vectors of bigint, center, and reverse order.
                let mut ct0i: Vec<BigInt> =
                    ct0_coeffs.iter().rev().map(|&x| BigInt::from(x)).collect();
                let mut ct1i: Vec<BigInt> =
                    ct1_coeffs.iter().rev().map(|&x| BigInt::from(x)).collect();
                let mut pk0i: Vec<BigInt> =
                    pk0_coeffs.iter().rev().map(|&x| BigInt::from(x)).collect();
                let mut pk1i: Vec<BigInt> =
                    pk1_coeffs.iter().rev().map(|&x| BigInt::from(x)).collect();

                let qi_bigint = BigInt::from(qi.modulus());

                reduce_and_center_coefficients_mut(&mut ct0i, &qi_bigint);
                reduce_and_center_coefficients_mut(&mut ct1i, &qi_bigint);
                reduce_and_center_coefficients_mut(&mut pk0i, &qi_bigint);
                reduce_and_center_coefficients_mut(&mut pk1i, &qi_bigint);

                // k0qi = -t^{-1} mod qi

                let koqi_u64 = qi.inv(qi.neg(t.modulus())).unwrap();
                let k0qi = BigInt::from(koqi_u64); // Do not need to center this

                // ki = k1 * k0qi
                let ki_poly = Polynomial::new(k1.clone()).scalar_mul(&k0qi);
                let ki = ki_poly.coefficients().to_vec();

                // Calculate ct0i_hat = pk0 * ui + e0i + ki
                let ct0i_hat = {
                    let pk0i_poly = Polynomial::new(pk0i.clone());
                    let u_poly = Polynomial::new(u.clone());
                    let pk0i_times_u = pk0i_poly.mul(&u_poly);
                    assert_eq!((pk0i_times_u.coefficients().len() as u64) - 1, 2 * (n - 1));

                    let e0_poly = Polynomial::new(e0.clone());
                    let ki_poly = Polynomial::new(ki.clone());
                    let e0_plus_ki = e0_poly.add(&ki_poly);
                    assert_eq!((e0_plus_ki.coefficients().len() as u64) - 1, n - 1);

                    pk0i_times_u.add(&e0_plus_ki).coefficients().to_vec()
                };
                assert_eq!((ct0i_hat.len() as u64) - 1, 2 * (n - 1));

                // Check whether ct0i_hat mod R_qi (the ring) is equal to ct0i
                let mut ct0i_hat_mod_rqi = ct0i_hat.clone();
                reduce_in_ring(&mut ct0i_hat_mod_rqi, &cyclo, &qi_bigint);
                assert_eq!(&ct0i, &ct0i_hat_mod_rqi);

                // Compute r2i numerator = ct0i - ct0i_hat and reduce/center the polynomial
                let ct0i_poly = Polynomial::new(ct0i.clone());
                let ct0i_hat_poly = Polynomial::new(ct0i_hat.clone());
                let ct0i_minus_ct0i_hat = ct0i_poly.sub(&ct0i_hat_poly).coefficients().to_vec();
                assert_eq!((ct0i_minus_ct0i_hat.len() as u64) - 1, 2 * (n - 1));
                let mut ct0i_minus_ct0i_hat_mod_zqi = ct0i_minus_ct0i_hat.clone();
                reduce_and_center_coefficients_mut(&mut ct0i_minus_ct0i_hat_mod_zqi, &qi_bigint);

                // Compute r2i as the quotient of numerator divided by the cyclotomic polynomial
                // to produce: (ct0i - ct0i_hat) / (x^N + 1) mod Z_qi. Remainder should be empty.
                let ct0i_minus_ct0i_hat_poly = Polynomial::new(ct0i_minus_ct0i_hat_mod_zqi.clone());
                let cyclo_poly = Polynomial::new(cyclo.clone());
                let (r2i_poly, r2i_rem_poly) = ct0i_minus_ct0i_hat_poly.div(&cyclo_poly).unwrap();
                let r2i = r2i_poly.coefficients().to_vec();
                let r2i_rem = r2i_rem_poly.coefficients().to_vec();
                assert!(r2i_rem.iter().all(|x| x.is_zero()));
                assert_eq!((r2i.len() as u64) - 1, n - 2); // Order(r2i) = N - 2

                // Assert that (ct0i - ct0i_hat) = (r2i * cyclo) mod Z_qi
                let r2i_poly = Polynomial::new(r2i.clone());
                let r2i_times_cyclo = r2i_poly.mul(&cyclo_poly).coefficients().to_vec();
                let mut r2i_times_cyclo_mod_zqi = r2i_times_cyclo.clone();
                reduce_and_center_coefficients_mut(&mut r2i_times_cyclo_mod_zqi, &qi_bigint);
                assert_eq!(&ct0i_minus_ct0i_hat_mod_zqi, &r2i_times_cyclo_mod_zqi);
                assert_eq!((r2i_times_cyclo.len() as u64) - 1, 2 * (n - 1));

                // Calculate r1i = (ct0i - ct0i_hat - r2i * cyclo) / qi mod Z_p. Remainder should be empty.
                let ct0i_minus_ct0i_hat_poly = Polynomial::new(ct0i_minus_ct0i_hat.clone());
                let r2i_times_cyclo_poly = Polynomial::new(r2i_times_cyclo.clone());
                let r1i_num = ct0i_minus_ct0i_hat_poly
                    .sub(&r2i_times_cyclo_poly)
                    .coefficients()
                    .to_vec();
                assert_eq!((r1i_num.len() as u64) - 1, 2 * (n - 1));

                let r1i_num_poly = Polynomial::new(r1i_num.clone());
                let qi_poly = Polynomial::new(vec![qi_bigint.clone()]);
                let (r1i_poly, r1i_rem_poly) = r1i_num_poly.div(&qi_poly).unwrap();
                let r1i = r1i_poly.coefficients().to_vec();
                let r1i_rem = r1i_rem_poly.coefficients().to_vec();
                assert!(r1i_rem.iter().all(|x| x.is_zero()));
                assert_eq!((r1i.len() as u64) - 1, 2 * (n - 1)); // Order(r1i) = 2*(N-1)
                let r1i_poly_check = Polynomial::new(r1i.clone());
                assert_eq!(
                    &r1i_num,
                    &r1i_poly_check.mul(&qi_poly).coefficients().to_vec()
                );

                // Assert that ct0i = ct0i_hat + r1i * qi + r2i * cyclo mod Z_p
                let r1i_poly = Polynomial::new(r1i.clone());
                let r1i_times_qi = r1i_poly.scalar_mul(&qi_bigint).coefficients().to_vec();
                let ct0i_hat_poly = Polynomial::new(ct0i_hat.clone());
                let r1i_times_qi_poly = Polynomial::new(r1i_times_qi.clone());
                let r2i_times_cyclo_poly = Polynomial::new(r2i_times_cyclo.clone());
                let mut ct0i_calculated = ct0i_hat_poly
                    .add(&r1i_times_qi_poly)
                    .add(&r2i_times_cyclo_poly)
                    .coefficients()
                    .to_vec();

                while ct0i_calculated.len() > 0 && ct0i_calculated[0].is_zero() {
                    ct0i_calculated.remove(0);
                }

                assert_eq!(&ct0i, &ct0i_calculated);

                // --------------------------------------------------- ct1i ---------------------------------------------------

                // Calculate ct1i_hat = pk1i * ui + e1i
                let ct1i_hat = {
                    let pk1i_poly = Polynomial::new(pk1i.clone());
                    let u_poly = Polynomial::new(u.clone());
                    let pk1i_times_u = pk1i_poly.mul(&u_poly);
                    assert_eq!((pk1i_times_u.coefficients().len() as u64) - 1, 2 * (n - 1));

                    let e1_poly = Polynomial::new(e1.clone());
                    pk1i_times_u.add(&e1_poly).coefficients().to_vec()
                };
                assert_eq!((ct1i_hat.len() as u64) - 1, 2 * (n - 1));

                // Check whether ct1i_hat mod R_qi (the ring) is equal to ct1i
                let mut ct1i_hat_mod_rqi = ct1i_hat.clone();
                reduce_in_ring(&mut ct1i_hat_mod_rqi, &cyclo, &qi_bigint);
                assert_eq!(&ct1i, &ct1i_hat_mod_rqi);

                // Compute p2i numerator = ct1i - ct1i_hat
                let ct1i_poly = Polynomial::new(ct1i.clone());
                let ct1i_hat_poly = Polynomial::new(ct1i_hat.clone());
                let ct1i_minus_ct1i_hat = ct1i_poly.sub(&ct1i_hat_poly).coefficients().to_vec();
                assert_eq!((ct1i_minus_ct1i_hat.len() as u64) - 1, 2 * (n - 1));
                let mut ct1i_minus_ct1i_hat_mod_zqi = ct1i_minus_ct1i_hat.clone();
                reduce_and_center_coefficients_mut(&mut ct1i_minus_ct1i_hat_mod_zqi, &qi_bigint);

                // Compute p2i as the quotient of numerator divided by the cyclotomic polynomial,
                // and reduce/center the resulting coefficients to produce:
                // (ct1i - ct1i_hat) / (x^N + 1) mod Z_qi. Remainder should be empty.
                let ct1i_minus_ct1i_hat_poly = Polynomial::new(ct1i_minus_ct1i_hat_mod_zqi.clone());
                let (p2i_poly, p2i_rem_poly) =
                    ct1i_minus_ct1i_hat_poly.div(&cyclo_poly.clone()).unwrap();
                let p2i = p2i_poly.coefficients().to_vec();
                let p2i_rem = p2i_rem_poly.coefficients().to_vec();
                assert!(p2i_rem.iter().all(|x| x.is_zero()));
                assert_eq!((p2i.len() as u64) - 1, n - 2); // Order(p2i) = N - 2

                // Assert that (ct1i - ct1i_hat) = (p2i * cyclo) mod Z_qi
                let p2i_poly = Polynomial::new(p2i.clone());
                let p2i_times_cyclo: Vec<BigInt> =
                    p2i_poly.mul(&cyclo_poly).coefficients().to_vec();
                let mut p2i_times_cyclo_mod_zqi = p2i_times_cyclo.clone();
                reduce_and_center_coefficients_mut(&mut p2i_times_cyclo_mod_zqi, &qi_bigint);
                assert_eq!(&ct1i_minus_ct1i_hat_mod_zqi, &p2i_times_cyclo_mod_zqi);
                assert_eq!((p2i_times_cyclo.len() as u64) - 1, 2 * (n - 1));

                // Calculate p1i = (ct1i - ct1i_hat - p2i * cyclo) / qi mod Z_p. Remainder should be empty.
                let ct1i_minus_ct1i_hat_poly = Polynomial::new(ct1i_minus_ct1i_hat.clone());
                let p2i_times_cyclo_poly = Polynomial::new(p2i_times_cyclo.clone());
                let p1i_num = ct1i_minus_ct1i_hat_poly
                    .sub(&p2i_times_cyclo_poly)
                    .coefficients()
                    .to_vec();
                assert_eq!((p1i_num.len() as u64) - 1, 2 * (n - 1));

                let p1i_num_poly = Polynomial::new(p1i_num.clone());
                let qi_poly = Polynomial::new(vec![BigInt::from(qi.modulus())]);
                let (p1i_poly, p1i_rem_poly) = p1i_num_poly.div(&qi_poly).unwrap();
                let p1i = p1i_poly.coefficients().to_vec();
                let p1i_rem = p1i_rem_poly.coefficients().to_vec();
                assert!(p1i_rem.iter().all(|x| x.is_zero()));
                assert_eq!((p1i.len() as u64) - 1, 2 * (n - 1)); // Order(p1i) = 2*(N-1)
                let p1i_poly_check = Polynomial::new(p1i.clone());
                assert_eq!(
                    &p1i_num,
                    &p1i_poly_check.mul(&qi_poly).coefficients().to_vec()
                );

                // Assert that ct1i = ct1i_hat + p1i * qi + p2i * cyclo mod Z_p
                let p1i_poly = Polynomial::new(p1i.clone());
                let p1i_times_qi = p1i_poly.scalar_mul(&qi_bigint).coefficients().to_vec();
                let ct1i_hat_poly = Polynomial::new(ct1i_hat.clone());
                let p1i_times_qi_poly = Polynomial::new(p1i_times_qi.clone());
                let p2i_times_cyclo_poly = Polynomial::new(p2i_times_cyclo.clone());
                let mut ct1i_calculated = ct1i_hat_poly
                    .add(&p1i_times_qi_poly)
                    .add(&p2i_times_cyclo_poly)
                    .coefficients()
                    .to_vec();

                while ct1i_calculated.len() > 0 && ct1i_calculated[0].is_zero() {
                    ct1i_calculated.remove(0);
                }

                assert_eq!(&ct1i, &ct1i_calculated);
                (i, r2i, r1i, k0qi, ct0i, ct1i, pk0i, pk1i, p1i, p2i)
            },
        )
        .collect();

        // Merge results into the `res` structure after parallel execution
        for (i, r2i, r1i, k0i, ct0i, ct1i, pk0i, pk1i, p1i, p2i) in results.into_iter() {
            res.r2is[i] = r2i;
            res.r1is[i] = r1i;
            res.k0is[i] = k0i;
            res.ct0is[i] = ct0i;
            res.ct1is[i] = ct1i;
            res.pk0is[i] = pk0i;
            res.pk1is[i] = pk1i;
            res.p1is[i] = p1i;
            res.p2is[i] = p2i;
        }

        // Set final result vectors
        res.u = u;
        res.e0 = e0;
        res.e1 = e1;
        res.k1 = k1;

        Ok(res)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fhe::bfv::{BfvParametersBuilder, Encoding, Plaintext, SecretKey};
    use num_bigint::BigInt;
    use rand::{rngs::StdRng, SeedableRng};
    use std::str::FromStr;

    fn setup_test_params() -> (Arc<BfvParameters>, SecretKey, PublicKey) {
        let params = BfvParametersBuilder::new()
            .set_degree(2048)
            .set_plaintext_modulus(1032193)
            .set_moduli(&[18014398492704769])
            .build_arc()
            .unwrap();

        let mut rng = StdRng::seed_from_u64(0); // Use deterministic seed
        let sk = SecretKey::random(&params, &mut rng);
        let pk = PublicKey::new(&sk, &mut rng);

        (params, sk, pk)
    }

    #[test]
    fn test_vector_lengths() {
        let vecs = InputValidationVectors::new(1, 2048);
        assert!(vecs.check_correct_lengths(1, 2048));
        assert!(!vecs.check_correct_lengths(2, 2048)); // Wrong moduli count
        assert!(!vecs.check_correct_lengths(1, 1024)); // Wrong degree
    }

    #[test]
    fn test_standard_form() {
        let vecs = InputValidationVectors::new(1, 2048);
        let p = BigInt::from_str(
            "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        )
        .unwrap();
        let std_form = vecs.standard_form(&p);

        // Check that all vectors are properly reduced
        assert!(std_form.u.iter().all(|x| x < &p));
        assert!(std_form.e0.iter().all(|x| x < &p));
        assert!(std_form.e1.iter().all(|x| x < &p));
        assert!(std_form.k1.iter().all(|x| x < &p));
    }

    #[test]
    fn test_vector_computation() {
        let (params, sk, pk) = setup_test_params();

        // Create a sample plaintext
        let mut message_data = vec![3u64; params.degree()];
        message_data[0] = 1;
        let pt = Plaintext::try_encode(&message_data, Encoding::poly(), &params).unwrap();

        // Use extended encryption to get the polynomial data
        let mut rng = StdRng::seed_from_u64(0);
        let (_ct, u_rns, e0_rns, e1_rns) = pk.try_encrypt_extended(&pt, &mut rng).unwrap();

        // Compute vectors
        let vecs =
            InputValidationVectors::compute(&pt, &u_rns, &e0_rns, &e1_rns, &_ct, &pk, &params)
                .unwrap();

        // Check dimensions
        assert!(vecs.check_correct_lengths(1, params.degree()));
    }

    #[test]
    fn test_vector_json_format() {
        let vecs = InputValidationVectors::new(1, 4); // Small size for testing
        let json = vecs.to_json();

        // Check all required fields are present
        let required_fields = [
            "pk0i", "pk1i", "u", "e0", "e1", "k1", "r2is", "r1is", "p2is", "p1is", "ct0is", "ct1is",
        ];

        for field in required_fields.iter() {
            assert!(json.get(field).is_some(), "Missing field: {}", field);
        }
    }
}
