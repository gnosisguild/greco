mod poly;

use fhe::bfv::{BfvParametersBuilder, Ciphertext, Encoding, Plaintext, PublicKey, SecretKey};
use fhe_math::{
    rq::{Context, Poly, Representation},
    zq::Modulus,
};
use fhe_traits::*;
use itertools::izip;
use num_bigint::BigInt;
use num_traits::{Num, Signed, ToPrimitive, Zero};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use serde_json::json;
use std::fs::File;
use std::io::Write;
use std::ops::Deref;
use std::path::Path;
use std::sync::Arc;
use std::vec;

use poly::*;

fn main() {
    // Set up the BFV parameters
    let N: u64 = 128;
    let plaintext_modulus: u64 = 65537;
    let moduli: Vec<u64> = vec![4503599625535489, 4503599626321921];

    let params = BfvParametersBuilder::new()
        .set_degree(N as usize)
        .set_plaintext_modulus(plaintext_modulus)
        .set_moduli(&moduli)
        .build_arc()
        .unwrap();

    // Extract plaintext modulus
    let t = Modulus::new(params.plaintext()).unwrap();

    // Use a seedable rng for experimental reproducibility
    let mut rng = StdRng::seed_from_u64(0);

    // Generate the secret and public keys
    let sk = SecretKey::random(&params, &mut rng);
    let pk = PublicKey::new(&sk, &mut rng);

    // Sample a message and encrypt
    // let m = t.random_vec(N as usize, &mut rng);
    let m: Vec<i64> = (-(N as i64 / 2)..(N as i64 / 2)).collect(); // m here is from lowest degree to largest as input into fhe.rs (REQUIRED)
    let pt = Plaintext::try_encode(&m, Encoding::poly(), &params).unwrap();
    let (ct, u_rns, e0_rns, e1_rns) = pk.try_encrypt_extended(&pt, &mut rng).unwrap();

    // Sanity check. m = Decrypt(ct)
    let m_decrypted = unsafe { t.center_vec_vt(&sk.try_decrypt(&ct).unwrap().value.into_vec()) };
    assert_eq!(m_decrypted, m);

    // Extract context
    let ctx = params.ctx_at_level(pt.level()).unwrap().clone();

    // Initialize zk proving modulus
    let p = BigInt::from_str_radix(
        "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        10,
    )
    .unwrap();

    let (r2is, r1is, k0is, ct0is, ct1is, pk0is, pk1is, p1is, p2is, u, e0, e1, k1) =
        compute_input_validation_vectors(&ctx, &t, &pt, &u_rns, &e0_rns, &e1_rns, &ct, &pk);

    // Create standard form versions with respect to p
    let (
        pk0is_std,
        pk1is_std,
        r2is_std,
        r1is_std,
        p2is_std,
        p1is_std,
        ct0is_std,
        ct1is_std,
        u_std,
        e0_std,
        e1_std,
        k1_std,
    ) = input_validation_vectors_standard_form(
        &pk0is, &pk1is, &r2is, &r1is, &p2is, &p1is, &ct0is, &ct1is, &u, &e0, &e1, &k1, &p,
    );

    // Create output json with standard form polynomials
    let json_data = json!({
        "pk0i": to_string_2d_vec(&pk0is_std),
        "pk1i": to_string_2d_vec(&pk1is_std),
        "u": to_string_1d_vec(&u_std),
        "e0": to_string_1d_vec(&e0_std),
        "e1": to_string_1d_vec(&e1_std),
        "k1": to_string_1d_vec(&k1_std),
        "r2is": to_string_2d_vec(&r2is_std),
        "r1is": to_string_2d_vec(&r1is_std),
        "p2is": to_string_2d_vec(&p2is_std),
        "p1is": to_string_2d_vec(&p1is_std),
        "ct0is": to_string_2d_vec(&ct0is_std),
        "ct1is": to_string_2d_vec(&ct1is_std)
    });

    let moduli_bigint: Vec<BigInt> = moduli.iter().map(|&qi| BigInt::from(qi)).collect();
    let moduli_bitsize = {
        if let Some(&max_value) = moduli.iter().max() {
            64 - max_value.leading_zeros()
        } else {
            0
        }
    };

    // Calculate bounds ---------------------------------------------------------------------

    // constraint. The coefficients of u, e0, e1 should be in the range [-⌈6σ⌋, ⌈6σ⌋]
    // where ⌈6σ⌋ is the upper bound of the discrete Gaussian distribution
    //
    // Note: the secret key in fhe.rs is sampled from a discrete gaussian distribution
    // rather than a ternary distribution as in bfv.py.
    let gauss_bound = BigInt::from(
        f64::ceil(6_f64 * f64::sqrt(params.variance() as f64))
            .to_i64()
            .unwrap(),
    );
    let u_bound = gauss_bound.clone();
    let e_bound = gauss_bound.clone();
    // Check centered bounds for u, e0, e1
    assert!(range_check_centered(&u, &-&u_bound, &u_bound));
    assert!(range_check_centered(&e0, &-&e_bound, &e_bound));
    assert!(range_check_centered(&e1, &-&e_bound, &e_bound));

    // Check assigned bounds for u, e0, e1
    assert!(range_check_standard(&u_std, &u_bound, &p));
    assert!(range_check_standard(&e0_std, &e_bound, &p));
    assert!(range_check_standard(&e1_std, &e_bound, &p));

    // constraint. The coefficients of k1 should be in the range [-(t-1)/2, (t-1)/2]
    let ptxt_bound = BigInt::from((t.modulus() - 1) / 2);
    let k1_bound = ptxt_bound.clone();
    assert!(range_check_centered(&k1, &-&k1_bound, &k1_bound));
    assert!(range_check_standard(&k1_std, &k1_bound, &p));

    // Calculate bounds and perform asserts for polynomials depending on each qi
    let mut pk_bounds: Vec<BigInt> = vec![BigInt::zero(); moduli.len()];
    let mut r2_bounds: Vec<BigInt> = vec![BigInt::zero(); moduli.len()];
    let mut r1_bounds: Vec<BigInt> = vec![BigInt::zero(); moduli.len()];
    let mut p2_bounds: Vec<BigInt> = vec![BigInt::zero(); moduli.len()];
    let mut p1_bounds: Vec<BigInt> = vec![BigInt::zero(); moduli.len()];
    for i in 0..moduli.len() {
        // constraint. The coefficients of ct0i and ct1i should be in the range [-(qi-1)/2, (qi-1)/2]
        let bound = (moduli_bigint[i].clone() - BigInt::from(1)) / BigInt::from(2);
        assert!(range_check_centered(&ct0is[i], &-&bound, &bound));
        assert!(range_check_centered(&ct1is[i], &-&bound, &bound));

        // constraint. The coefficients of pk0i and pk1i should be in range [-(qi-1)/2 , (qi-1)/2]
        pk_bounds[i] = bound.clone();
        assert!(range_check_centered(
            &pk0is[i],
            &-&pk_bounds[i],
            &pk_bounds[i]
        ));
        assert!(range_check_centered(
            &pk1is[i],
            &-&pk_bounds[i],
            &pk_bounds[i]
        ));
        assert!(range_check_standard(&pk0is_std[i], &pk_bounds[i], &p));
        assert!(range_check_standard(&pk1is_std[i], &pk_bounds[i], &p));

        // constraint. The coefficients of r2i should be in the range [-(qi-1)/2, (qi-1)/2]
        r2_bounds[i] = bound.clone();
        assert!(range_check_centered(
            &r2is[i],
            &-&r2_bounds[i],
            &r2_bounds[i]
        ));
        assert!(range_check_standard(&r2is_std[i], &r2_bounds[i], &p));

        // constraint. The coefficients of (ct0i - ct0i_hat - r2i * cyclo) / qi = r1i should be in the range
        // $[
        //      \frac{- ((N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|)}{q_i},
        //      \frac{   (N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}| }{q_i}
        // ]$
        let r1i_bound =
            (BigInt::from(&N + 2) * &bound + &gauss_bound + &ptxt_bound * BigInt::abs(&k0is[i]))
                / &moduli_bigint[i];
        r1_bounds[i] = r1i_bound.clone();
        assert!(range_check_centered(
            &r1is[i],
            &-&r1_bounds[i],
            &r1_bounds[i]
        ));
        assert!(range_check_standard(&r1is_std[i], &r1i_bound, &p));

        // constraint. The coefficients of p2 should be in the range [-(qi-1)/2, (qi-1)/2]
        p2_bounds[i] = bound.clone();
        assert!(range_check_centered(
            &p2is[i],
            &-&p2_bounds[i],
            &p2_bounds[i]
        ));
        assert!(range_check_standard(&p2is_std[i], &p2_bounds[i], &p));

        // constraint. The coefficients of (ct0i - ct0i_hat - p2i * cyclo) / qi = p1i should be in the range
        // $[
        //      \frac{- ((N+2) \cdot \frac{q_i - 1}{2} + B)}{q_i},
        //      \frac{   (N+2) \cdot \frac{q_i - 1}{2} + B }{q_i}
        // ]$
        let p1i_bound = (BigInt::from(&N + 2) * &bound + &gauss_bound) / &moduli_bigint[i];
        p1_bounds[i] = p1i_bound.clone();
        assert!(range_check_centered(
            &p1is[i],
            &-&p1_bounds[i],
            &p1_bounds[i]
        ));
        assert!(range_check_standard(&p1is_std[i], &p1i_bound, &p));
    }

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
    let zeroes_json = create_zeroes_json(params.degree(), moduli.len());
    write_json_to_file(&output_path, &filename_zeroes, &zeroes_json);

    let filename_constants = format!(
        "pk_enc_constants_{}_{}x{}_{}.rs",
        N,
        moduli.len(),
        moduli_bitsize,
        t.modulus()
    );

    write_constants_to_file(
        &N,
        &pk_bounds,
        &e_bound,
        &u_bound,
        &r1_bounds,
        &r2_bounds,
        &p1_bounds,
        &p2_bounds,
        &k1_bound,
        &moduli_bigint,
        &k0is,
        &filename_constants,
    )
}

/// Create the centered validation vectors necessary for creating an input validation proof according to Greco.
/// For more information, please see https://eprint.iacr.org/2024/594.
///
/// # Arguments
///
/// * `ctx` - Context object from fhe.rs holding information about elements in Rq.
/// * `t` - Plaintext modulus object.
/// * `pt` - Plaintext from fhe.rs.
/// * `u_rns` - Private polynomial used in ciphertext sampled from secret key distribution.
/// * `e0_rns` - Error polynomial used in ciphertext sampled from error distribution.
/// * `e1_rns` - Error polynomioal used in cihpertext sampled from error distribution.
/// * `ct` - Ciphertext from fhe.rs.
/// * `pk` - Public Key from fhe.re.
///
fn compute_input_validation_vectors(
    ctx: &Arc<Context>,
    t: &Modulus,
    pt: &Plaintext,
    u_rns: &Poly,
    e0_rns: &Poly,
    e1_rns: &Poly,
    ct: &Ciphertext,
    pk: &PublicKey,
) -> (
    Vec<Vec<BigInt>>,
    Vec<Vec<BigInt>>,
    Vec<BigInt>,
    Vec<Vec<BigInt>>,
    Vec<Vec<BigInt>>,
    Vec<Vec<BigInt>>,
    Vec<Vec<BigInt>>,
    Vec<Vec<BigInt>>,
    Vec<Vec<BigInt>>,
    Vec<BigInt>,
    Vec<BigInt>,
    Vec<BigInt>,
    Vec<BigInt>,
) {
    let N: u64 = ctx.degree as u64;

    // Calculate k1 (independent of qi), center and reverse
    let q_mod_t = (ctx.modulus() % t.modulus()).to_u64().unwrap(); // [q]_t
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
            .center_vec_vt(u_rns_copy.coefficients().row(0).as_slice().unwrap())
            .iter()
            .map(|&x| BigInt::from(x))
            .rev()
            .collect()
    };

    let e0: Vec<BigInt> = unsafe {
        ctx.moduli_operators()[0]
            .center_vec_vt(e0_rns_copy.coefficients().row(0).as_slice().unwrap())
            .iter()
            .map(|&x| BigInt::from(x))
            .rev()
            .collect()
    };

    let e1: Vec<BigInt> = unsafe {
        ctx.moduli_operators()[0]
            .center_vec_vt(e1_rns_copy.coefficients().row(0).as_slice().unwrap())
            .iter()
            .map(|&x| BigInt::from(x))
            .rev()
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
    let mut r2is: Vec<Vec<BigInt>> = vec![Vec::new(); num_moduli];
    let mut r1is: Vec<Vec<BigInt>> = vec![Vec::new(); num_moduli];
    let mut k0is: Vec<BigInt> = vec![BigInt::zero(); num_moduli];
    let mut ct0is: Vec<Vec<BigInt>> = vec![Vec::new(); num_moduli];
    let mut ct0is_hat: Vec<Vec<BigInt>> = vec![Vec::new(); num_moduli];
    let mut ct1is: Vec<Vec<BigInt>> = vec![Vec::new(); num_moduli];
    let mut ct1is_hat: Vec<Vec<BigInt>> = vec![Vec::new(); num_moduli];
    let mut pk0is: Vec<Vec<BigInt>> = vec![Vec::new(); num_moduli];
    let mut pk1is: Vec<Vec<BigInt>> = vec![Vec::new(); num_moduli];
    let mut p1is: Vec<Vec<BigInt>> = vec![Vec::new(); num_moduli];
    let mut p2is: Vec<Vec<BigInt>> = vec![Vec::new(); num_moduli];

    // Initialize iterators for results calculation
    let moduli_operators = ctx.moduli_operators();
    let ct0_iter = ct0.coefficients();
    let ct1_iter = ct1.coefficients();
    let pk0_iter = pk0.coefficients();
    let pk1_iter = pk1.coefficients();
    let zipped: Vec<_> = izip!(
        moduli_operators,
        ct0_iter.rows(),
        ct1_iter.rows(),
        pk0_iter.rows(),
        pk1_iter.rows()
    )
    .collect();

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
        Vec<BigInt>,
        Vec<BigInt>,
    )> = zipped
        .into_par_iter()
        .enumerate()
        .map(
            |(i, (qi, ct0_coeffs, ct1_coeffs, pk0_coeffs, pk1_coeffs))| {
                // --------------------------------------------------- ct0i ---------------------------------------------------

                // Convert to vectors of bigint, center, and reverse order.
                let mut ct0i: Vec<BigInt> =
                    ct0_coeffs.iter().map(|&x| BigInt::from(x)).rev().collect();
                let mut ct1i: Vec<BigInt> =
                    ct1_coeffs.iter().map(|&x| BigInt::from(x)).rev().collect();
                let mut pk0i: Vec<BigInt> =
                    pk0_coeffs.iter().map(|&x| BigInt::from(x)).rev().collect();
                let mut pk1i: Vec<BigInt> =
                    pk1_coeffs.iter().map(|&x| BigInt::from(x)).rev().collect();

                let qi_bigint = BigInt::from(qi.modulus());

                reduce_and_center_coefficients_mut(&mut ct0i, &qi_bigint);
                reduce_and_center_coefficients_mut(&mut ct1i, &qi_bigint);
                reduce_and_center_coefficients_mut(&mut pk0i, &qi_bigint);
                reduce_and_center_coefficients_mut(&mut pk1i, &qi_bigint);

                // k0qi = -t^{-1} mod qi
                let koqi_u64 = qi.inv(qi.neg(t.modulus())).unwrap();
                let k0qi = BigInt::from(koqi_u64); // Do not need to center this

                // ki = k1 * k0qi
                let ki = poly_scalar_mul(&k1, &k0qi);

                // Calculate ct0i_hat = pk0 * ui + e0i + ki
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

                // --------------------------------------------------- ct1i ---------------------------------------------------

                // Calculate ct1i_hat = pk1i * ui + e1i
                let ct1i_hat = {
                    let pk1i_times_u = poly_mul(&pk1i, &u);
                    assert_eq!((pk1i_times_u.len() as u64) - 1, 2 * (N - 1));

                    poly_add(&pk1i_times_u, &e1)
                };
                assert_eq!((ct1i_hat.len() as u64) - 1, 2 * (N - 1));

                // Check whether ct1i_hat mod R_qi (the ring) is equal to ct1i
                let mut ct1i_hat_mod_rqi = ct1i_hat.clone();
                reduce_in_ring(&mut ct1i_hat_mod_rqi, &cyclo, &qi_bigint);
                assert_eq!(&ct1i, &ct1i_hat_mod_rqi);

                // Compute p2i numerator = ct1i - ct1i_hat
                let ct1i_minus_ct1i_hat = poly_sub(&ct1i, &ct1i_hat);
                assert_eq!((ct1i_minus_ct1i_hat.len() as u64) - 1, 2 * (N - 1));
                let mut ct1i_minus_ct1i_hat_mod_zqi = ct1i_minus_ct1i_hat.clone();
                reduce_and_center_coefficients_mut(&mut ct1i_minus_ct1i_hat_mod_zqi, &qi_bigint);

                // Compute p2i as the quotient of numerator divided by the cyclotomic polynomial,
                // and reduce/center the resulting coefficients to produce:
                // (ct1i - ct1i_hat) / (x^N + 1) mod Z_qi. Remainder should be empty.
                let (p2i, p2i_rem) = poly_div(&ct1i_minus_ct1i_hat_mod_zqi, &cyclo.clone());
                assert!(p2i_rem.is_empty());
                assert_eq!((p2i.len() as u64) - 1, N - 2); // Order(p2i) = N - 2

                // Assert that (ct1i - ct1i_hat) = (p2i * cyclo) mod Z_qi
                let p2i_times_cyclo: Vec<BigInt> = poly_mul(&p2i, &cyclo);
                let mut p2i_times_cyclo_mod_zqi = p2i_times_cyclo.clone();
                reduce_and_center_coefficients_mut(&mut p2i_times_cyclo_mod_zqi, &qi_bigint);
                assert_eq!(&ct1i_minus_ct1i_hat_mod_zqi, &p2i_times_cyclo_mod_zqi);
                assert_eq!((p2i_times_cyclo.len() as u64) - 1, 2 * (N - 1));

                // Calculate p1i = (ct1i - ct1i_hat - p2i * cyclo) / qi mod Z_p. Remainder should be empty.
                let p1i_num = poly_sub(&ct1i_minus_ct1i_hat, &p2i_times_cyclo);
                assert_eq!((p1i_num.len() as u64) - 1, 2 * (N - 1));

                let (p1i, p1i_rem) = poly_div(&p1i_num, &[BigInt::from(qi.modulus())]);
                assert!(p1i_rem.is_empty());
                assert_eq!((p1i.len() as u64) - 1, 2 * (N - 1)); // Order(p1i) = 2*(N-1)
                assert_eq!(&p1i_num, &poly_mul(&p1i, &[qi_bigint.clone()]));

                // Assert that ct1i = ct1i_hat + p1i * qi + p2i * cyclo mod Z_p
                let p1i_times_qi = poly_scalar_mul(&p1i, &qi_bigint);
                let mut ct1i_calculated =
                    poly_add(&poly_add(&ct1i_hat, &p1i_times_qi), &p2i_times_cyclo);

                while ct1i_calculated.len() > 0 && ct1i_calculated[0].is_zero() {
                    ct1i_calculated.remove(0);
                }

                assert_eq!(&ct1i, &ct1i_calculated);

                /*
                println!("qi = {:?}\n", &qi_bigint);
                println!("ct0i = {:?}\n", &ct0i);
                println!("k0qi = {:?}\n", &k0qi);
                println!("pk0 = Polynomial({:?})\n", &pk0i);
                println!("pk1 = Polynomial({:?})\n", &pk1i);
                println!("ki = {:?}\n", &ki);
                println!("ct0i_hat_mod_rqi = {:?}\n", &ct0i_hat_mod_rqi);
                */

                (
                    i, r2i, r1i, k0qi, ct0i, ct0i_hat, ct1i, ct1i_hat, pk0i, pk1i, p1i, p2i,
                )
            },
        )
        .collect();

    // println!("Completed creation of polynomials!");

    // Aggregate results into global vectors
    for (i, r2i, r1i, k0i, ct0i, ct0i_hat, ct1i, ct1i_hat, pk0i, pk1i, p1i, p2i) in
        results.into_iter()
    {
        r2is[i] = r2i;
        r1is[i] = r1i;
        k0is[i] = k0i;
        ct0is[i] = ct0i;
        ct0is_hat[i] = ct0i_hat;
        ct1is[i] = ct1i;
        ct1is_hat[i] = ct1i_hat;
        pk0is[i] = pk0i;
        pk1is[i] = pk1i;
        p1is[i] = p1i;
        p2is[i] = p2i;
    }

    (
        r2is, r1is, k0is, ct0is, ct1is, pk0is, pk1is, p1is, p2is, u, e0, e1, k1,
    )
}

/// Assign and return all of the centered input validation vectors to the ZKP modulus `p`.
///
/// # Arguments
///
/// * `pk0is` - Centered coefficients of first public key object for each RNS modulus
/// * `pk1is` - Centered coefficients of second public key object for each RNS modulus
/// * `r2is` - Centered coefficients of r2 for each RNS modulus
/// * `r1is` - Centered coefficients of r1 for each RNS modulus
/// * `p2is` - Centered coefficients of p2 for each RNS modulus
/// * `p1is` - Centered coefficients of p1 for each RNS modulus
/// * `ct0is` - Centered coefficients of first ciphertext object for each RNS modulus
/// * `ct1is` - Centered coefficients of second ciphertext object for each RNS modulus
/// * `u` - Centered coefficients of secret polynomial used during encryption (sampled from secret key distribution)
/// * `e0` - Centered coefficients of error polynomial used during encryption (sampled from error distribution)
/// * `e1` - Centered coefficients of error polynomial used during encryption (sampled from error distribution)
/// * `k1` - Centered coefficients of [Q*m] mod t
/// * `p` - ZKP modulus
///
pub fn input_validation_vectors_standard_form(
    pk0is: &[Vec<BigInt>],
    pk1is: &[Vec<BigInt>],
    r2is: &[Vec<BigInt>],
    r1is: &[Vec<BigInt>],
    p2is: &[Vec<BigInt>],
    p1is: &[Vec<BigInt>],
    ct0is: &[Vec<BigInt>],
    ct1is: &[Vec<BigInt>],
    u: &[BigInt],
    e0: &[BigInt],
    e1: &[BigInt],
    k1: &[BigInt],
    p: &BigInt,
) -> (
    Vec<Vec<BigInt>>,
    Vec<Vec<BigInt>>,
    Vec<Vec<BigInt>>,
    Vec<Vec<BigInt>>,
    Vec<Vec<BigInt>>,
    Vec<Vec<BigInt>>,
    Vec<Vec<BigInt>>,
    Vec<Vec<BigInt>>,
    Vec<BigInt>,
    Vec<BigInt>,
    Vec<BigInt>,
    Vec<BigInt>,
) {
    (
        reduce_coefficients_2d(pk0is, p),
        reduce_coefficients_2d(pk1is, p),
        reduce_coefficients_2d(r2is, p),
        reduce_coefficients_2d(r1is, p),
        reduce_coefficients_2d(p2is, p),
        reduce_coefficients_2d(p1is, p),
        reduce_coefficients_2d(ct0is, p),
        reduce_coefficients_2d(ct1is, p),
        reduce_coefficients(u, p),
        reduce_coefficients(e0, p),
        reduce_coefficients(e1, p),
        reduce_coefficients(k1, p),
    )
}

fn write_constants_to_file(
    n: &u64,
    pk_bounds: &[BigInt],
    e_bound: &BigInt,
    u_bound: &BigInt,
    r1_bounds: &[BigInt],
    r2_bounds: &[BigInt],
    p1_bounds: &[BigInt],
    p2_bounds: &[BigInt],
    k1_bound: &BigInt,
    qi_constants: &[BigInt],
    k0i_constants: &[BigInt],
    output_file: &str,
) {
    // Set the output file path
    let output_path = Path::new("src")
        .join("constants")
        .join("pk_enc_constants")
        .join(output_file);

    let mut file = File::create(output_path).expect("Unable to create file");

    // Writing the constants to the file
    writeln!(file, "/// `N` is the degree of the cyclotomic polynomial defining the ring `Rq = Zq[X]/(X^N + 1)`.")
        .expect("Unable to write to file");
    writeln!(file, "pub const N: usize = {};", n).expect("Unable to write to file");

    let pk_bound_str = pk_bounds
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<String>>()
        .join(", ");
    writeln!(file, "/// The coefficients of the polynomial `pk0is` and `pk1is` should exist in the interval `[-PK_BOUND, PK_BOUND]`.")
        .expect("Unable to write to file");
    writeln!(
        file,
        "pub const PK_BOUND: [u64; {}] = [{}];",
        pk_bounds.len(),
        pk_bound_str
    )
    .expect("Unable to write to file");

    writeln!(file, "/// The coefficients of the polynomial `pk1is` should exist in the interval `[-PK0_BOUND, PK0_BOUND]`.")
        .expect("Unable to write to file");

    writeln!(file, "/// The coefficients of the polynomial `e` should exist in the interval `[-E_BOUND, E_BOUND]` where `E_BOUND` is the upper bound of the gaussian distribution with 𝜎 = 3.2.")
        .expect("Unable to write to file");
    writeln!(file, "pub const E_BOUND: u64 = {};", e_bound).expect("Unable to write to file");

    writeln!(file, "/// The coefficients of the polynomial `s` should exist in the interval `[-S_BOUND, S_BOUND]`.")
        .expect("Unable to write to file");
    writeln!(file, "pub const U_BOUND: u64 = {};", u_bound).expect("Unable to write to file");

    let r1_bounds_str = r1_bounds
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<String>>()
        .join(", ");
    writeln!(file, "/// The coefficients of the polynomials `r1is` should exist in the interval `[-R1_BOUND[i], R1_BOUND[i]]` where `R1_BOUND[i]` is equal to `(qi-1)/2`.")
        .expect("Unable to write to file");
    writeln!(
        file,
        "pub const R1_BOUNDS: [u64; {}] = [{}];",
        r1_bounds.len(),
        r1_bounds_str
    )
    .expect("Unable to write to file");

    let r2_bounds_str = r2_bounds
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<String>>()
        .join(", ");
    writeln!(file, "/// The coefficients of the polynomials `r2is` should exist in the interval `[-R2_BOUND[i], R2_BOUND[i]]` where `R2_BOUND[i]` is equal to $\\frac{{(N+2) \\cdot \\frac{{q_i - 1}}{{2}} + B + \\frac{{t - 1}}{{2}} \\cdot |K_{{0,i}}|}}{{q_i}}`.")
        .expect("Unable to write to file");
    writeln!(
        file,
        "pub const R2_BOUNDS: [u64; {}] = [{}];",
        r2_bounds.len(),
        r2_bounds_str
    )
    .expect("Unable to write to file");

    let p1_bounds_str = p1_bounds
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<String>>()
        .join(", ");
    writeln!(file, "/// The coefficients of the polynomials `p1is` should exist in the interval `[-P1_BOUND[i], P1_BOUND[i]]` where `P1_BOUND[i]` is equal to (((qis[i] - 1) / 2) * (n + 2) + b ) / qis[i].")
        .expect("Unable to write to file");
    writeln!(
        file,
        "pub const P1_BOUNDS: [u64; {}] = [{}];",
        p1_bounds.len(),
        p1_bounds_str
    )
    .expect("Unable to write to file");

    let p2_bounds_str = p2_bounds
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<String>>()
        .join(", ");
    writeln!(file, "/// The coefficients of the polynomials `p2is` should exist in the interval `[-P2_BOUND[i], P2_BOUND[i]]` where `P2_BOUND[i]` is equal to (qis[i] - 1) / 2.")
        .expect("Unable to write to file");
    writeln!(
        file,
        "pub const P2_BOUNDS: [u64; {}] = [{}];",
        p2_bounds.len(),
        p2_bounds_str
    )
    .expect("Unable to write to file");

    writeln!(file, "/// The coefficients of `k1` should exist in the interval `[-K1_BOUND, K1_BOUND]` where `K1_BOUND` is equal to `(t-1)/2`.")
        .expect("Unable to write to file");
    writeln!(file, "pub const K1_BOUND: u64 = {};", k1_bound).expect("Unable to write to file");

    let qis_str = qi_constants
        .iter()
        .map(|x| format!("\"{}\"", x))
        .collect::<Vec<String>>()
        .join(", ");
    writeln!(file, "/// List of scalars `qis` such that `qis[i]` is the modulus of the i-th CRT basis of `q` (ciphertext space modulus).")
        .expect("Unable to write to file");
    writeln!(
        file,
        "pub const QIS: [&str; {}] = [{}];",
        qi_constants.len(),
        qis_str
    )
    .expect("Unable to write to file");

    let k0is_str = k0i_constants
        .iter()
        .map(|x| format!("\"{}\"", x))
        .collect::<Vec<String>>()
        .join(", ");
    writeln!(file, "/// List of scalars `k0is` such that `k0i[i]` is equal to the negative of the multiplicative inverses of t mod qi.")
        .expect("Unable to write to file");
    writeln!(
        file,
        "pub const K0IS: [&str; {}] = [{}];",
        k0i_constants.len(),
        k0is_str
    )
    .expect("Unable to write to file");
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

fn create_zeroes_json(degree: usize, moduli_len: usize) -> serde_json::Value {
    json!({
        "pk0i": vec![vec![String::from("0"); degree]; moduli_len],
        "pk1i": vec![vec![String::from("0"); degree]; moduli_len],
        "u": vec![String::from("0"); degree],
        "e0": vec![String::from("0"); degree],
        "e1": vec![String::from("0"); degree],
        "k1": vec![String::from("0"); degree],
        "r2is": vec![vec![String::from("0"); degree]; moduli_len],
        "r1is": vec![vec![String::from("0"); degree]; moduli_len],
        "p2is": vec![vec![String::from("0"); degree]; moduli_len],
        "p1is": vec![vec![String::from("0"); degree]; moduli_len],
        "ct0is": vec![vec![String::from("0"); degree]; moduli_len],
        "ct1is": vec![vec![String::from("0"); degree]; moduli_len]
    })
}
