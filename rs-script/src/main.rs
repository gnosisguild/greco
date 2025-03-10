mod poly;

use fhe::bfv::{
    BfvParameters, BfvParametersBuilder, Ciphertext, Encoding, Plaintext, PublicKey, SecretKey,
};
use fhe_math::{
    rq::{Poly, Representation},
    zq::Modulus,
};
use fhe_traits::*;
use itertools::izip;

use num_bigint::BigInt;
use num_traits::{Num, Signed, ToPrimitive, Zero};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rayon::iter::{ParallelBridge, ParallelIterator};
use serde_json::json;
use std::fs::File;
use std::io::Write;
use std::ops::Deref;
use std::path::Path;
use std::sync::Arc;
use std::vec;

use poly::*;

/// Set of vectors for input validation of a ciphertext
#[derive(Clone, Debug)]
pub struct KeyGenValidationVectors {
    pk0is: Vec<Vec<BigInt>>,
    pk1is: Vec<Vec<BigInt>>,
    r1is: Vec<Vec<BigInt>>,
    r2is: Vec<Vec<BigInt>>,
    s: Vec<BigInt>,
    e: Vec<BigInt>,
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
    /// TODO Set the right degrees for these polynomials 
    pub fn new(num_moduli: usize, degree: usize) -> Self {
        InputValidationVectors {
            pk0is: vec![vec![BigInt::zero(); degree]; num_moduli],
            pk1is: vec![vec![BigInt::zero(); degree]; num_moduli],
            r1is: vec![vec![BigInt::zero(); 2 * (degree - 1) + 1]; num_moduli],
            r2is: vec![vec![BigInt::zero(); degree - 1]; num_moduli],
            s: vec![BigInt::zero(); degree],
            e: vec![BigInt::zero(); degree],
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
            r1is: reduce_coefficients_2d(&self.r1is, p),
            r2is: reduce_coefficients_2d(&self.r2is, p),
            s: reduce_coefficients(&self.s, p),
            e: reduce_coefficients(&self.e, p),
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
            "s": to_string_1d_vec(&self.s),
            "e": to_string_1d_vec(&self.e),
            "r2is": to_string_2d_vec(&self.r2is),
            "r1is": to_string_2d_vec(&self.r1is),
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
            check_2d_lengths(&self.r1is, num_moduli, 2 * (degree - 1)),
            check_2d_lengths(&self.r2is, num_moduli, degree - 2),
            // 1D vector checks
            check_1d_lengths(&self.s, degree),
            check_1d_lengths(&self.e, degree),
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
        s_rns: &Poly,
        e_rns: &Poly,
        pk: &PublicKey,
    ) -> Result<InputValidationVectors, Box<dyn std::error::Error>> {
        // Get context, plaintext modulus, and degree
        let params = &pk.par;
        let ctx = params.ctx_at_level(pt.level())?;
        let t = Modulus::new(params.plaintext())?;
        let N: u64 = ctx.degree as u64;

        let q_mod_t = (ctx.modulus() % t.modulus())
            .to_u64()
            .ok_or_else(|| "Cannot convert BigInt to u64.".to_string())?; // [q]_t

        // Extract single vectors of u, e1, and e2 as Vec<BigInt>, center and reverse
        let mut s_rns_copy = s_rns.clone();
        let mut e_rns_copy = e_rns.clone();

        s_rns_copy.change_representation(Representation::PowerBasis);
        e_rns_copy.change_representation(Representation::PowerBasis);


        let s: Vec<BigInt> = unsafe {
            ctx.moduli_operators()[0]
                .center_vec_vt(
                    s_rns_copy
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

        let e: Vec<BigInt> = unsafe {
            ctx.moduli_operators()[0]
                .center_vec_vt(
                    e_rns_copy
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

        let mut pk0: Poly = pk.c.c[0].clone();
        let mut pk1: Poly = pk.c.c[1].clone();
        pk0.change_representation(Representation::PowerBasis);
        pk1.change_representation(Representation::PowerBasis);

        // Create cyclotomic polynomial x^N + 1

        let mut cyclo = vec![BigInt::from(0u64); (N + 1) as usize];

        cyclo[0] = BigInt::from(1u64); // x^N term
        cyclo[N as usize] = BigInt::from(1u64); // x^0 term

        // Print
        /*
        println!("m = {:?}\n", &m);
        println!("k1 = {:?}\n", &k1);
        println!("u = {:?}\n", &u);
        println!("e0 = {:?}\n", &e0);
        println!("e1 = {:?}\n", &e1);
         */

        // Initialize matrices to store results
        let num_moduli = ctx.moduli().len();
        let mut res = InputValidationVectors::new(num_moduli, N as usize);

        // Perform the main computation logic
        // TODO to correct this..
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
            pk0.coefficients().rows(),
            pk1.coefficients().rows()
        )
        .enumerate()
        .par_bridge()
        .map(
            |(i, (qi, pk0_coeffs, pk1_coeffs))| {
                // --------------------------------------------------- ct0i ---------------------------------------------------

                // Convert to vectors of bigint, center, and reverse order.
                let mut pk0i: Vec<BigInt> =
                    pk0_coeffs.iter().rev().map(|&x| BigInt::from(x)).collect();
                let mut pk1i: Vec<BigInt> =
                    pk1_coeffs.iter().rev().map(|&x| BigInt::from(x)).collect();

                let qi_bigint = BigInt::from(qi.modulus());

                reduce_and_center_coefficients_mut(&mut pk0i, &qi_bigint);
                reduce_and_center_coefficients_mut(&mut pk1i, &qi_bigint);

                // TODO To correct the caculations happening below..
                // Calculate ct0i_hat = -pk1 * s + e0i 
                let ct0i_hat = {
                    let pk0i_times_u = poly_mul(&pk0i, &u);
                    assert_eq!((pk0i_times_u.len() as u64) - 1, 2 * (N - 1));
                    let e0_plus_ki = poly_add(&e0, &ki);
                    assert_eq!((e0_plus_ki.len() as u64) - 1, N - 1);

                    poly_add(&pk0i_times_u, &e0_plus_ki)
                };
                assert_eq!((ct0i_hat.len() as u64) - 1, 2 * (N - 1));

                // Check whether ct0i_hat mod R_qi (the ring) is equal to ct0i
                let mut ct0i_hat_mod_rqi = ct0i_hat.clone();
                reduce_in_ring(&mut ct0i_hat_mod_rqi, &cyclo, &qi_bigint);
                assert_eq!(&ct0i, &ct0i_hat_mod_rqi);

                // Compute r2i numerator = ct0i - ct0i_hat and reduce/center the polynomial
                let ct0i_minus_ct0i_hat = poly_sub(&ct0i, &ct0i_hat);
                assert_eq!((ct0i_minus_ct0i_hat.len() as u64) - 1, 2 * (N - 1));
                let mut ct0i_minus_ct0i_hat_mod_zqi = ct0i_minus_ct0i_hat.clone();
                reduce_and_center_coefficients_mut(&mut ct0i_minus_ct0i_hat_mod_zqi, &qi_bigint);

                // Compute r2i as the quotient of numerator divided by the cyclotomic polynomial
                // to produce: (ct0i - ct0i_hat) / (x^N + 1) mod Z_qi. Remainder should be empty.
                let (r2i, r2i_rem) = poly_div(&ct0i_minus_ct0i_hat_mod_zqi, &cyclo);
                assert!(r2i_rem.is_empty());
                assert_eq!((r2i.len() as u64) - 1, N - 2); // Order(r2i) = N - 2

                // Assert that (ct0i - ct0i_hat) = (r2i * cyclo) mod Z_qi
                let r2i_times_cyclo = poly_mul(&r2i, &cyclo);
                let mut r2i_times_cyclo_mod_zqi = r2i_times_cyclo.clone();
                reduce_and_center_coefficients_mut(&mut r2i_times_cyclo_mod_zqi, &qi_bigint);
                assert_eq!(&ct0i_minus_ct0i_hat_mod_zqi, &r2i_times_cyclo_mod_zqi);
                assert_eq!((r2i_times_cyclo.len() as u64) - 1, 2 * (N - 1));

                // Calculate r1i = (ct0i - ct0i_hat - r2i * cyclo) / qi mod Z_p. Remainder should be empty.
                let r1i_num = poly_sub(&ct0i_minus_ct0i_hat, &r2i_times_cyclo);
                assert_eq!((r1i_num.len() as u64) - 1, 2 * (N - 1));

                let (r1i, r1i_rem) = poly_div(&r1i_num, &[qi_bigint.clone()]);
                assert!(r1i_rem.is_empty());
                assert_eq!((r1i.len() as u64) - 1, 2 * (N - 1)); // Order(r1i) = 2*(N-1)
                assert_eq!(&r1i_num, &poly_mul(&r1i, &[qi_bigint.clone()]));

                // Assert that ct0i = ct0i_hat + r1i * qi + r2i * cyclo mod Z_p
                let r1i_times_qi = poly_scalar_mul(&r1i, &qi_bigint);
                let mut ct0i_calculated =
                    poly_add(&poly_add(&ct0i_hat, &r1i_times_qi), &r2i_times_cyclo);

                while ct0i_calculated.len() > 0 && ct0i_calculated[0].is_zero() {
                    ct0i_calculated.remove(0);
                }

                assert_eq!(&ct0i, &ct0i_calculated);

                (i, r2i, r1i, pk0i, pk1i)
            },
        )
        .collect();

        // println!("Completed creation of polynomials!");

        // Merge results into the `res` structure after parallel execution
        for (i, r2i, r1i, k0i, ct0i, ct1i, pk0i, pk1i, p1i, p2i) in results.into_iter() {
            res.r2is[i] = r2i;
            res.r1is[i] = r1i;
            res.pk0is[i] = pk0i;
            res.pk1is[i] = pk1i;
        }

        // Set final result vectors
        res.u = u;
        res.e = e;
        Ok(res)
    }
}

/// The `InputValidationBounds` struct holds the bounds for various vectors and polynomials used in the input validation process.
/// These bounds are calculated from a set of BFV encryption parameters and represent limits on the values of different fields
/// to ensure that the inputs remain within valid ranges during operations.
/// TODO: To Check why we had pk here
#[derive(Clone, Debug)]
pub struct InputValidationBounds {
    s: BigInt,
    e: BigInt,
    pk: Vec<BigInt>,
    r1_low: Vec<BigInt>,
    r1_up: Vec<BigInt>,
    r2: Vec<BigInt>,
}

impl InputValidationBounds {
    /// Checks the constraints of the input validation vectors against the bounds stored in `InputValidationBounds`.
    ///
    /// # Arguments
    ///
    /// * `vecs` - A reference to `InputValidationVectors`, which contains the vectors to be validated.
    /// * `p` - The prime modulus used in the encryption scheme.
    ///
    /// This function checks whether the coefficients of the vectors `u`, `e0`, `e1`, `k1`, and others are within
    /// the specified ranges, using both centered and standard range checks. It asserts that the vectors stay within
    /// these predefined bounds.
    pub fn check_constraints(&self, vecs: &InputValidationVectors, p: &BigInt) {
        let vecs_std = vecs.standard_form(p);

        // constraint. The coefficients of u, e0, e1 should be in the range [-⌈6σ⌋, ⌈6σ⌋]
        // where ⌈6σ⌋ is the upper bound of the discrete Gaussian distribution
        // TODO check why we had both centered and standard checks..
        assert!(range_check_centered(&vecs.s, &-&self.s, &self.s));
        assert!(range_check_centered(&vecs.e, &-&self.e, &self.e));
        assert!(range_check_standard(&vecs_std.s, &self.s, &p));
        assert!(range_check_standard(&vecs_std.e, &self.e, &p));


        // Perform asserts for polynomials depending on each qi
        for i in 0..self.r2.len() {
            // constraint. The coefficients of pk0i and pk1i should be in range [-(qi-1)/2 , (qi-1)/2]
            assert!(range_check_centered(
                &vecs.pk0is[i],
                &-&self.pk[i],
                &self.pk[i]
            ));
            assert!(range_check_centered(
                &vecs.pk1is[i],
                &-&self.pk[i],
                &self.pk[i]
            ));
            assert!(range_check_standard(&vecs_std.pk0is[i], &self.pk[i], &p));
            assert!(range_check_standard(&vecs_std.pk1is[i], &self.pk[i], &p));

            // constraint. The coefficients of r2i should be in the range [-(qi-1)/2, (qi-1)/2]
            assert!(range_check_centered(
                &vecs.r2is[i],
                &-&self.r2[i],
                &self.r2[i]
            ));
            assert!(range_check_standard(&vecs_std.r2is[i], &self.r2[i], &p));

            // constraint. The coefficients of (ct0i - ct0i_hat - r2i * cyclo) / qi = r1i should be in the range
            // $[
            //      \frac{ \frac{-(t - 1)}{2} \cdot |K_{0,i}| - ((N \cdot B +2) \cdot \frac{q_i - 1}{2} + B )}{q_i},
            //      \frac{ \frac{t - 1}{2} \cdot |K_{0,i}| +  (N \cdot B+2) \cdot \frac{q_i - 1}{2} + B }{q_i}
            // ]$
            assert!(range_check_centered(
                &vecs.r1is[i],
                &self.r1_low[i],
                &self.r1_up[i]
            ));
            assert!(range_check_standard_2bounds(
                &vecs_std.r1is[i],
                &self.r1_low[i],
                &self.r1_up[i],
                &p
            ));

    }

    /// Compute the input validation bounds from a set of BFV encryption parameters.
    ///
    /// # Arguments
    ///
    /// * `params` - A reference to the BFV parameters.
    /// * `level` - The encryption level, which determines the number of moduli used.
    ///
    /// # Returns
    ///
    /// A new `InputValidationBounds` instance containing the bounds for vectors and polynomials
    /// based on the BFV parameters and the specified level.
    pub fn compute(
        params: &Arc<BfvParameters>,
        level: usize,
    ) -> Result<InputValidationBounds, Box<dyn std::error::Error>> {
        // Get cyclotomic degree and context at provided level
        let N = BigInt::from(params.degree());
        let t = BigInt::from(params.plaintext());
        let ctx = params.ctx_at_level(level)?;

        // Note: the secret key in fhe.rs is sampled from a discrete gaussian distribution
        // rather than a ternary distribution as in bfv.py.
        let gauss_bound = BigInt::from(
            f64::ceil(6_f64 * f64::sqrt(params.variance() as f64))
                .to_i64()
                .ok_or_else(|| "Failed to convert variance to i64".to_string())?,
        );
        let s_bound = gauss_bound.clone();
        let e_bound = gauss_bound.clone();

        //Note we have two different variables for lower bound and upper bound, as in the case
        //where the plaintext modulus is even, the lower bound cannot be calculated by just
        //negating the upper bound. For instance, if t = 8, then the lower bound will be -4 and the
        //upper bound will be 3
        //
        // Calculate qi-based bounds
        // TODO correct these bounds
        let num_moduli = ctx.moduli().len();
        let mut pk_bounds: Vec<BigInt> = vec![BigInt::zero(); num_moduli];
        let mut r2_bounds: Vec<BigInt> = vec![BigInt::zero(); num_moduli];
        let mut r1_low_bounds: Vec<BigInt> = vec![BigInt::zero(); num_moduli];
        let mut r1_up_bounds: Vec<BigInt> = vec![BigInt::zero(); num_moduli];
        for (i, qi) in ctx.moduli_operators().iter().enumerate() {
            let qi_bigint = BigInt::from(qi.modulus());
            let qi_bound = (&qi_bigint - BigInt::from(1)) / BigInt::from(2);

            // Calculate the k0qi for the bounds (these are also constant wrt BFV params)
            let k0qi = BigInt::from(
                qi.inv(qi.neg(params.plaintext()))
                    .ok_or_else(|| "Failed to calculate modulus inverse for k0qi".to_string())?,
            );

            pk_bounds[i] = qi_bound.clone();
            r2_bounds[i] = qi_bound.clone();

            r1_low_bounds[i] = (&ptxt_low_bound * BigInt::abs(&k0qi)
                - &((&N * &gauss_bound + 2) * &qi_bound + &gauss_bound))
                / &qi_bigint;
            r1_up_bounds[i] = (&ptxt_up_bound * BigInt::abs(&k0qi)
                + ((&N * &gauss_bound + 2) * &qi_bound + &gauss_bound))
                / &qi_bigint;

        }

        Ok(InputValidationBounds {
            s: s_bound,
            e: e_bound,
            pk: pk_bounds,
            r1_low: r1_low_bounds,
            r1_up: r1_up_bounds,
            r2: r2_bounds,
        })
    }

    /// Writes the input validation bounds to a file that can be imported as a Rust module.
    ///
    /// # Arguments
    ///
    /// * `params` - Reference to BFV parameters to extract context information.
    /// * `output_file` - The path where the output constants should be saved.
    ///
    /// This function calculates certain constants like `k0i` values for each modulus `qi` and writes the bounds and other
    /// relevant constants in a Rust-friendly format to the file specified by `output_file`.
    fn to_file(
        &self,
        params: &Arc<BfvParameters>,
        output_file: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let level = params.moduli().len() - self.r2.len();
        let ctx = params.ctx_at_level(level)?;

        // Set the output file path
        let output_path = Path::new("src")
            .join("constants")
            .join("pk_enc_constants")
            .join(output_file);

        let mut file = File::create(output_path)?;

        // Writing the constants to the file
        writeln!(file, "/// `N` is the degree of the cyclotomic polynomial defining the ring `Rq = Zq[X]/(X^N + 1)`.")?;
        writeln!(file, "pub const N: usize = {};", params.degree())?;

        let pk_bound_str = self
            .pk
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join(", ");
        writeln!(file, "/// The coefficients of the polynomial `pk0is` and `pk1is` should exist in the interval `[-PK_BOUND, PK_BOUND]`.")?;
        writeln!(
            file,
            "pub const PK_BOUND: [u64; {}] = [{}];",
            self.pk.len(),
            pk_bound_str
        )?;

        writeln!(file, "/// The coefficients of the polynomial `e` should exist in the interval `[-E_BOUND, E_BOUND]` where `E_BOUND` is the upper bound of the gaussian distribution with 𝜎 = 3.2.")?;
        writeln!(file, "pub const E_BOUND: u64 = {};", self.e)?;

        writeln!(file, "/// The coefficients of the polynomial `s` should exist in the interval `[-S_BOUND, S_BOUND]`.")?;
        writeln!(file, "pub const S_BOUND: u64 = {};", self.s)?;

        let r1_low_bounds_str = self
            .r1_low
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join(", ");

        let r1_up_bounds_str = self
            .r1_up
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join(", ");
        writeln!(file, "/// The coefficients of the polynomials `r1is` should exist in the interval `[R1_LOW_BOUNDS[i], R1_UP_BOUNDS[i]]` where R1_LOW_BOUNDS is equal to $\\frac{{\\frac{{-(t - 1)}}{{2}} \\cdot |K_{{0,i}}| - (N \\cdot B +2) \\cdot \\frac{{q_i - 1}}{{2}} + B}}{{q_i}}` and `R1_UP_BOUNDS[i]` is equal to `$\\frac{{\\frac{{(t - 1)}}{{2}} \\cdot |K_{{0,i}}| + (N \\cdot +2) \\cdot \\frac{{q_i - 1}}{{2}} + B}}{{q_i}}` .")?;
        writeln!(
            file,
            "pub const R1_LOW_BOUNDS: [i64; {}] = [{}];",
            self.r1_low.len(),
            r1_low_bounds_str
        )?;
        writeln!(
            file,
            "pub const R1_UP_BOUNDS: [u64; {}] = [{}];",
            self.r1_up.len(),
            r1_up_bounds_str
        )?;

        let r2_bounds_str = self
            .r2
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join(", ");
        writeln!(file, "/// The coefficients of the polynomials `r2is` should exist in the interval `[-R2_BOUND[i], R2_BOUND[i]]` where `R2_BOUND[i]` is equal to `(qi-1)/2`.")?;
        writeln!(
            file,
            "pub const R2_BOUNDS: [u64; {}] = [{}];",
            self.r2.len(),
            r2_bounds_str
        )?;


        let qis_str = ctx
            .moduli()
            .iter()
            .map(|x| format!("\"{}\"", x))
            .collect::<Vec<String>>()
            .join(", ");
        writeln!(file, "/// List of scalars `qis` such that `qis[i]` is the modulus of the i-th CRT basis of `q` (ciphertext space modulus).")?;
        writeln!(
            file,
            "pub const QIS: [&str; {}] = [{}];",
            ctx.moduli().len(),
            qis_str
        )?;

        Ok(())
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Set up the BFV parameters

    let N: u64 = 1024;

    //let plaintext_modulus: u64 = 65537;
    let plaintext_modulus: u64 = 2048;

    // let moduli: Vec<u64> = vec![
    //   1152921504606584833,
    //   1152921504598720513,
    //   1152921504597016577,
    //   1152921504595968001,
    //   1152921504595640321,
    //   1152921504593412097,
    //   1152921504592822273,
    //   1152921504592429057,
    //   1152921504589938689,
    //   1152921504586530817,
    //   1152921504585547777,
    //   1152921504583647233,
    //   1152921504581877761,
    //   1152921504581419009,
    //   1152921504580894721
    // ];

    let moduli: Vec<u64> = vec![4503599625535489, 4503599626321921];
    //let moduli: Vec<u64> = vec![4503599625535489, 4503599626321921];
    //let moduli: Vec<u64> = vec![1038337,18014398492704769,4503599625535489, 4503599626321921];
    //let moduli: Vec<u64> = vec![1038337];
    //let moduli: Vec<u64> = vec![
    //18014398492704769
    //];

    //let moduli_sizes = [20];

    let params = BfvParametersBuilder::new()
        .set_degree(N as usize)
        .set_plaintext_modulus(plaintext_modulus)
        .set_moduli(&moduli)
        //.set_moduli_sizes(&moduli_sizes)
        .build_arc()?;

    // Extract plaintext modulus
    let t = Modulus::new(params.plaintext())?;

    // Use a seedable rng for experimental reproducibility
    let mut rng = StdRng::seed_from_u64(0);

    // Generate the secret and public keys
    let sk = SecretKey::random(&params, &mut rng);

    let pk = PublicKey::new(&sk, &mut rng);

    //Sample a message and encrypt
    let m = t.random_vec(N as usize, &mut rng);
    //let m: Vec<i64> = (-(N as i64 / 2)..(N as i64 / 2)).collect(); // m here is from lowest degree to largest as input into fhe.rs (REQUIRED)
    let pt = Plaintext::try_encode(&m, Encoding::poly(), &params)?;

    // Extract context
    let ctx = params.ctx_at_level(pt.level())?.clone();
    let (ct, u_rns, e0_rns, e1_rns) = pk.try_encrypt_extended(&pt, &mut rng)?;

    // Sanity check. m = Decrypt(ct)

    //let m_decrypted = unsafe { t.center_vec_vt(&sk.try_decrypt(&ct)?.value.into_vec()) };
    let m_decrypted = sk.try_decrypt(&ct)?.value.into_vec();
    assert_eq!(m_decrypted, m);

    let moduli = ctx.moduli();
    // Initialize zk proving modulus
    let p = BigInt::from_str_radix(
        "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        10,
    )?;

    // Compute input validation vectors
    let res = InputValidationVectors::compute(&pt, &u_rns, &e0_rns, &e1_rns, &ct, &pk)?;

    // Create output json with standard form polynomials
    let json_data = res.standard_form(&p).to_json();

    // Calculate bounds ---------------------------------------------------------------------
    let bounds = InputValidationBounds::compute(&params, pt.level())?;

    // Check the constraints
    bounds.check_constraints(&res, &p);

    let moduli_bitsize = {
        if let Some(&max_value) = ctx.moduli().iter().max() {
            64 - max_value.leading_zeros()
        } else {
            0
        }
    };

    // Write out files ----------------------------------------------------------------------
    let output_path = Path::new("src").join("data").join("pk_enc_data");

    // Generate filename and write file
    let filename = format!(
        "pk_enc_{}_{}x{}_{}.json",
        N,
        moduli.len(),
        moduli_bitsize,
        t.modulus()
    );
    write_json_to_file(&output_path, &filename, &json_data);

    // Generate zeros filename and write file
    let filename_zeroes = format!(
        "pk_enc_{}_{}x{}_{}_zeroes.json",
        N,
        moduli.len(),
        moduli_bitsize,
        t.modulus()
    );
    let zeroes_json = InputValidationVectors::new(moduli.len(), params.degree()).to_json();
    write_json_to_file(&output_path, &filename_zeroes, &zeroes_json);

    let filename_constants = format!(
        "pk_enc_constants_{}_{}x{}_{}.rs",
        N,
        moduli.len(),
        moduli_bitsize,
        t.modulus()
    );
    bounds.to_file(&params, &filename_constants)?;

    Ok(())
}

fn to_string_1d_vec(poly: &Vec<BigInt>) -> Vec<String> {
    poly.iter().map(|coef| coef.to_string()).collect()
}

fn to_string_2d_vec(poly: &Vec<Vec<BigInt>>) -> Vec<Vec<String>> {
    poly.iter().map(|row| to_string_1d_vec(row)).collect()
}

/// Writes the given JSON data to a file in the specified output path.
///
/// # Arguments
///
/// * `output_path` - A reference to the base directory path where the file will be created.
/// * `filename` - The name of the file to create.
/// * `json_data` - A reference to the JSON data to be written into the file.
///
/// # Panics
///
/// This function will panic if the file cannot be created or if writing to the file fails.
fn write_json_to_file(output_path: &Path, filename: &str, json_data: &serde_json::Value) {
    let file_path = output_path.join(filename);
    let mut file = File::create(file_path).expect("Unable to create file");
    file.write_all(serde_json::to_string_pretty(json_data).unwrap().as_bytes())
        .expect("Unable to write data");
}
