use fhe::bfv::{BfvParametersBuilder, Encoding, Plaintext, PublicKey, SecretKey};
use fhe_math::rq::{Poly, Representation};
use fhe_math::zq::Modulus;
use fhe_traits::*;
use itertools::izip;
use num_bigint::BigInt;
use num_traits::{Num, ToPrimitive, Zero};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rayon::prelude::*;
use serde_json::json;
use std::fs::File;
use std::io::Write;
use std::ops::Deref;
use std::path::Path;
use std::vec;

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
    let m: Vec<i64> = (-64..64).collect(); // m here is from lowest degree to largest as input into fhe.rs (REQUIRED)
    let pt = Plaintext::try_encode(&m, Encoding::poly(), &params).unwrap();
    let (ct, mut u_rns, mut e0_rns, mut e1_rns) = pk.try_encrypt_extended(&pt, &mut rng).unwrap();

    // Sanity check. m = Decrypt(ct)
    let m_decrypted = unsafe { t.center_vec_vt(&sk.try_decrypt(&ct).unwrap().value.into_vec()) };
    assert_eq!(m_decrypted, m);

    // Extract context
    let ctx = params.ctx_at_level(pt.level()).unwrap().clone();

    // Calculate k1 (independent of qi), center and reverse
    let q_mod_t = (ctx.modulus() % t.modulus()).to_u64().unwrap(); // [q]_t
    let mut k1_u64 = pt.value.deref().to_vec(); // m
    t.scalar_mul_vec(&mut k1_u64, q_mod_t); // k1 = [q*m]_t
    let mut k1: Vec<BigInt> = k1_u64.iter().map(|&x| BigInt::from(x)).rev().collect();
    reduce_and_center_coefficients_mut(&mut k1, &BigInt::from(t.modulus()));

    // Extract single vectors of u, e1, and e2 as Vec<BigInt>, center and reverse
    u_rns.change_representation(Representation::PowerBasis);
    e0_rns.change_representation(Representation::PowerBasis);
    e1_rns.change_representation(Representation::PowerBasis);
    let u: Vec<BigInt> = unsafe {
        ctx.moduli_operators()[0]
            .center_vec_vt(u_rns.coefficients().row(0).as_slice().unwrap())
            .iter()
            .map(|&x| BigInt::from(x))
            .rev()
            .collect()
    };

    let e0: Vec<BigInt> = unsafe {
        ctx.moduli_operators()[0]
            .center_vec_vt(e0_rns.coefficients().row(0).as_slice().unwrap())
            .iter()
            .map(|&x| BigInt::from(x))
            .rev()
            .collect()
    };

    let e1: Vec<BigInt> = unsafe {
        ctx.moduli_operators()[0]
            .center_vec_vt(e1_rns.coefficients().row(0).as_slice().unwrap())
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

    // Initialize zk proving modulus
    let p = BigInt::from_str_radix(
        "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        10,
    )
    .unwrap();

    // Print
    /*
    println!("m = {:?}\n", &m);
    println!("k1 = {:?}\n", &k1);
    println!("u = {:?}\n", &u);
    println!("e0 = {:?}\n", &e0);
    println!("e1 = {:?}\n", &e1);
     */

    // Initialize matrices to store results
    let num_moduli = params.moduli().len();
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

    println!("Completed creation of polynomials!");

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

    // Create standard form versions with respect to p
    let pk0is_assigned = reduce_coefficients_2d(&pk0is, &p);
    let pk1is_assigned = reduce_coefficients_2d(&pk1is, &p);
    let u_assigned = reduce_coefficients(&u, &p);
    let e0_assigned = reduce_coefficients(&e0, &p);
    let e1_assigned = reduce_coefficients(&e1, &p);
    let k1_assigned = reduce_coefficients(&k1, &p);
    let r2is_assigned = reduce_coefficients_2d(&r2is, &p);
    let r1is_assigned = reduce_coefficients_2d(&r1is, &p);
    let p2is_assigned = reduce_coefficients_2d(&p2is, &p);
    let p1is_assigned = reduce_coefficients_2d(&p1is, &p);
    let ct0is_assigned = reduce_coefficients_2d(&ct0is, &p);
    let ct1is_assigned = reduce_coefficients_2d(&ct1is, &p);

    // Create output json with standard form polynomials
    let json_data = json!({
        "pk0i": to_string_2d_vec(&pk0is_assigned),
        "pk1i": to_string_2d_vec(&pk1is_assigned),
        "u": to_string_1d_vec(&u_assigned),
        "e0": to_string_1d_vec(&e0_assigned),
        "e1": to_string_1d_vec(&e1_assigned),
        "k1": to_string_1d_vec(&k1_assigned),
        "r2is": to_string_2d_vec(&r2is_assigned),
        "r1is": to_string_2d_vec(&r1is_assigned),
        "p2is": to_string_2d_vec(&p2is_assigned),
        "p1is": to_string_2d_vec(&p1is_assigned),
        "ct0is": to_string_2d_vec(&ct0is_assigned),
        "ct1is": to_string_2d_vec(&ct1is_assigned)
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

    let key_bound = BigInt::from(
        f64::ceil(6_f64 * f64::sqrt(params.variance() as f64))
            .to_i64()
            .unwrap(),
    ); // round(6*sigma)
    let ptxt_bound = BigInt::from((t.modulus() - 1) / 2);
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
        assert!(range_check_standard(&pk0is_assigned[i], &pk_bounds[i], &p));
        assert!(range_check_standard(&pk1is_assigned[i], &pk_bounds[i], &p));

        // constraint. The coefficients of r2i should be in the range [-(qi-1)/2, (qi-1)/2]
        r2_bounds[i] = bound.clone();
        assert!(range_check_centered(
            &r2is[i],
            &-&r2_bounds[i],
            &r2_bounds[i]
        ));
        assert!(range_check_standard(&r2is_assigned[i], &r2_bounds[i], &p));

        // constraint. The coefficients of (ct0i - ct0i_hat - r2i * cyclo) / qi = r1i should be in the range
        // $[
        //      \frac{- ((N+2) \cdot \frac{q_i - 1}{2} + B +\frac{t - 1}{2} \cdot |K_{0,i}|)}{q_i},
        //      \frac{(N+2) \cdot \frac{q_i - 1}{2} + B + \frac{t - 1}{2} \cdot |K_{0,i}|}{q_i}
        // ]$
        let r1i_bound = (BigInt::from(&N + 2) * &bound + &key_bound + &ptxt_bound * &k0is[i])
            / &moduli_bigint[i];
        r1_bounds[i] = r1i_bound.clone();
        assert!(range_check_centered(
            &r1is[i],
            &-&r1_bounds[i],
            &r1_bounds[i]
        ));
        assert!(range_check_standard(&r1is_assigned[i], &r1i_bound, &p));

        // constraint. The coefficients of p2 should be in the range [-(qi-1)/2, (qi-1)/2]
        p2_bounds[i] = bound.clone();
        assert!(range_check_centered(
            &p2is[i],
            &-&p2_bounds[i],
            &p2_bounds[i]
        ));
        assert!(range_check_standard(&p2is_assigned[i], &p2_bounds[i], &p));
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

    // write_constants_to_file(
    //     N,
    //     pk_bound,
    //     b,
    //     r1_bounds,
    //     r2_bounds,
    //     p1_bounds,
    //     p2_bounds,
    //     k1_bound,
    //     qi_constants,
    //     k0i_constants,
    //     &filename_constants,
    // )
}

fn write_constants_to_file(
    n: usize,
    pk_bound: &[BigInt],
    b: BigInt,
    r1_bounds: &[BigInt],
    r2_bounds: &[BigInt],
    p1_bounds: &[BigInt],
    p2_bounds: &[BigInt],
    k1_bound: BigInt,
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

    let pk_bound_str = pk_bound
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<String>>()
        .join(", ");
    writeln!(file, "/// The coefficients of the polynomial `pk0is` and `pk1is` should exist in the interval `[-PK_BOUND, PK_BOUND]`.")
        .expect("Unable to write to file");
    writeln!(
        file,
        "pub const PK_BOUND: [u64; {}] = [{}];",
        pk_bound.len(),
        pk_bound_str
    )
    .expect("Unable to write to file");

    writeln!(file, "/// The coefficients of the polynomial `pk1is` should exist in the interval `[-PK0_BOUND, PK0_BOUND]`.")
        .expect("Unable to write to file");

    writeln!(file, "/// The coefficients of the polynomial `e` should exist in the interval `[-E_BOUND, E_BOUND]` where `E_BOUND` is the upper bound of the gaussian distribution with 𝜎 = 3.2.")
        .expect("Unable to write to file");
    writeln!(file, "pub const E_BOUND: u64 = {};", b).expect("Unable to write to file");

    writeln!(file, "/// The coefficients of the polynomial `s` should exist in the interval `[-S_BOUND, S_BOUND]`.")
        .expect("Unable to write to file");
    writeln!(file, "pub const U_BOUND: u64 = {};", 1).expect("Unable to write to file");

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

/// Adds two polynomials represented as vectors of `BigInt` coefficients in descending order of powers.
///
/// This function aligns two polynomials of potentially different lengths and adds their coefficients.
/// It assumes that polynomials are represented from leading degree to degree zero, even if the
/// coefficient at degree zero is zero. Leading zeros are not removed to keep the order of the
/// polynomial correct, which in Greco's case is necessary so that the order can be checked.
///
/// # Arguments
///
/// * `poly1` - Coefficients of the first polynomial, from highest to lowest degree.
/// * `poly2` - Coefficients of the second polynomial, from highest to lowest degree.
///
/// # Returns
///
/// A vector of `BigInt` coefficients representing the sum of the two polynomials.
fn poly_add(poly1: &[BigInt], poly2: &[BigInt]) -> Vec<BigInt> {
    // Determine the new length and create extended polynomials
    let max_length = std::cmp::max(poly1.len(), poly2.len());
    let mut extended_poly1 = vec![BigInt::zero(); max_length];
    let mut extended_poly2 = vec![BigInt::zero(); max_length];

    // Copy original coefficients into extended vectors
    extended_poly1[max_length - poly1.len()..].clone_from_slice(poly1);
    extended_poly2[max_length - poly2.len()..].clone_from_slice(poly2);

    // Add the coefficients
    let mut result = vec![BigInt::zero(); max_length];
    for i in 0..max_length {
        result[i] = &extended_poly1[i] + &extended_poly2[i];
    }

    result
}

/// Negates the coefficients of a polynomial represented as a slice of `BigInt` coefficients.
///
/// This function creates a new polynomial where each coefficient is the negation of the corresponding
/// coefficient in the input polynomial.
///
/// # Arguments
///
/// * `poly` - A slice of `BigInt` representing the coefficients of the polynomial, with the highest
///   degree term at index 0 and the constant term at the end.
///
/// # Returns
///
/// A vector of `BigInt` representing the polynomial with negated coefficients, with the same degree
/// order as the input polynomial.
fn poly_neg(poly: &[BigInt]) -> Vec<BigInt> {
    poly.iter().map(|x| -x).collect()
}

/// Subtracts one polynomial from another, both represented as slices of `BigInt` coefficients in descending order.
///
/// This function subtracts the second polynomial (`poly2`) from the first polynomial (`poly1`). It does so
/// by first negating the coefficients of `poly2` and then adding the result to `poly1`.
///
/// # Arguments
///
/// * `poly1` - A slice of `BigInt` representing the coefficients of the first polynomial (minuend), with
///   the highest degree term at index 0 and the constant term at the end.
/// * `poly2` - A slice of `BigInt` representing the coefficients of the second polynomial (subtrahend), with
///   the highest degree term at index 0 and the constant term at the end.
///
/// # Returns
///
/// A vector of `BigInt` representing the coefficients of the resulting polynomial after subtraction.
fn poly_sub(poly1: &[BigInt], poly2: &[BigInt]) -> Vec<BigInt> {
    poly_add(poly1, &poly_neg(poly2))
}

/// Multiplies two polynomials represented as slices of `BigInt` coefficients naively.
///
/// Given two polynomials `poly1` and `poly2`, where each polynomial is represented by a slice of
/// coefficients, this function computes their product. The order of coefficients (ascending or
/// descending powers) should be the same for both input polynomials. The resulting polynomial is
/// returned as a vector of `BigInt` coefficients in the same order as the inputs.
///
/// # Arguments
///
/// * `poly1` - A slice of `BigInt` representing the coefficients of the first polynomial.
/// * `poly2` - A slice of `BigInt` representing the coefficients of the second polynomial.
///
/// # Returns
///
/// A vector of `BigInt` representing the coefficients of the resulting polynomial after multiplication,
/// in the same order as the input polynomials.
fn poly_mul(poly1: &[BigInt], poly2: &[BigInt]) -> Vec<BigInt> {
    let product_len = poly1.len() + poly2.len() - 1;
    let mut product = vec![BigInt::zero(); product_len];

    for i in 0..poly1.len() {
        for j in 0..poly2.len() {
            product[i + j] += &poly1[i] * &poly2[j];
        }
    }

    product
}

/// Divides one polynomial by another, returning the quotient and remainder, with both polynomials
/// represented by vectors of `BigInt` coefficients in descending order of powers.
///
/// Given two polynomials `dividend` and `divisor`, where each polynomial is represented by a vector
/// of coefficients in descending order of powers (i.e., the coefficient at index `i` corresponds
/// to the term of degree `n - i`, where `n` is the degree of the polynomial), this function computes
/// their quotient and remainder. The quotient and remainder are also represented in descending order
/// of powers.
///
/// # Arguments
///
/// * `dividend` - A slice of `BigInt` representing the coefficients of the dividend polynomial.
/// * `divisor` - A slice of `BigInt` representing the coefficients of the divisor polynomial. The leading
///   coefficient (highest degree term) must be non-zero.
///
/// # Returns
///
/// A tuple containing two vectors of `BigInt`:
/// * The first vector represents the quotient polynomial, with coefficients in descending order of powers.
/// * The second vector represents the remainder polynomial, also in descending order of powers.
///
/// # Panics
///
/// This function will panic if the divisor is empty or if the leading coefficient of the divisor is zero.
fn poly_div(dividend: &[BigInt], divisor: &[BigInt]) -> (Vec<BigInt>, Vec<BigInt>) {
    assert!(
        !divisor.is_empty() && !divisor[0].is_zero(),
        "Leading coefficient of divisor cannot be zero"
    );

    let mut quotient = vec![BigInt::zero(); dividend.len() - divisor.len() + 1];
    let mut remainder = dividend.to_vec();

    for i in 0..quotient.len() {
        let coeff = &remainder[i] / &divisor[0];
        quotient[i] = coeff.clone();

        for j in 0..divisor.len() {
            remainder[i + j] = &remainder[i + j] - &divisor[j] * &coeff;
        }
    }

    while remainder.len() > 0 && remainder[0].is_zero() {
        remainder.remove(0);
    }

    (quotient, remainder)
}

/// Multiplies each coefficient of a polynomial by a scalar.
///
/// This function takes a polynomial represented as a vector of `BigInt` coefficients and multiplies each
/// coefficient by a given scalar.
///
/// # Arguments
///
/// * `poly` - A slice of `BigInt` representing the coefficients of the polynomial, with the highest degree term
///   at index 0 and the constant term at the end.
/// * `scalar` - A `BigInt` representing the scalar by which each coefficient of the polynomial will be multiplied.
///
/// # Returns
///
/// A vector of `BigInt` representing the polynomial with each coefficient multiplied by the scalar, maintaining
/// the same order of coefficients as the input polynomial.
fn poly_scalar_mul(poly: &[BigInt], scalar: &BigInt) -> Vec<BigInt> {
    poly.iter().map(|coeff| coeff * scalar).collect()
}

/// Reduces the coefficients of a polynomial by dividing it with a cyclotomic polynomial
/// and updating the coefficients with the remainder.
///
/// This function performs a polynomial long division of the input polynomial (represented by
/// `coefficients`) by the given cyclotomic polynomial (represented by `cyclo`). It replaces
/// the original coefficients with the coefficients of the remainder from this division.
///
/// # Arguments
///
/// * `coefficients` - A mutable reference to a `Vec<BigInt>` containing the coefficients of
///   the polynomial to be reduced. The coefficients are in descending order of degree,
///   i.e., the first element is the coefficient of the highest degree term.
/// * `cyclo` - A slice of `BigInt` representing the coefficients of the cyclotomic polynomial.
///   The coefficients are also in descending order of degree.
///
/// # Panics
///
/// This function will panic if the remainder length exceeds the degree of the cyclotomic polynomial,
/// which would indicate an issue with the division operation.
fn reduce_coefficients_by_cyclo(coefficients: &mut Vec<BigInt>, cyclo: &[BigInt]) {
    // Perform polynomial long division, assuming poly_div returns (quotient, remainder)
    let (_, remainder) = poly_div(&coefficients, cyclo);

    let N = cyclo.len() - 1;
    let mut out: Vec<BigInt> = vec![BigInt::zero(); N];

    // Calculate the starting index in `out` where the remainder should be copied
    let start_idx = N - remainder.len();

    // Copy the remainder into the `out` vector starting from `start_idx`
    out[start_idx..].clone_from_slice(&remainder);

    // Resize the original `coefficients` vector to fit the result and copy the values
    coefficients.clear();
    coefficients.extend(out);
}

/// Reduces a number modulo a prime modulus and centers it.
///
/// This function takes an arbitrary number and reduces it modulo the specified prime modulus.
/// After reduction, the number is adjusted to be within the symmetric range
/// [−(modulus−1)/2, (modulus−1)/2]. If the number is already within this range, it remains unchanged.
///
/// # Parameters
///
/// - `x`: A reference to a `BigInt` representing the number to be reduced and centered.
/// - `modulus`: A reference to the prime modulus `BigInt` used for reduction.
/// - `half_modulus`: A reference to a `BigInt` representing half of the modulus used to center the coefficient.
///
/// # Returns
///
/// - A `BigInt` representing the reduced and centered number.
fn reduce_and_center(x: &BigInt, modulus: &BigInt, half_modulus: &BigInt) -> BigInt {
    // Calculate the remainder ensuring it's non-negative
    let mut r = x % modulus;
    if r < BigInt::zero() {
        r += modulus;
    }

    // Adjust the remainder if it is greater than half_modulus
    if r > *half_modulus {
        r -= modulus;
    }

    r
}

/// Reduces and centers polynomial coefficients modulo a prime modulus.
///
/// This function iterates over a mutable slice of polynomial coefficients, reducing each coefficient
/// modulo a given prime modulus and adjusting the result to be within the symmetric range
/// [−(modulus−1)/2, (modulus−1)/2].
///
/// # Parameters
///
/// - `coefficients`: A mutable slice of `BigInt` coefficients to be reduced and centered.
/// - `modulus`: A prime modulus `BigInt` used for reduction and centering.
///
/// # Panics
///
/// - Panics if `modulus` is zero due to division by zero.
fn reduce_and_center_coefficients_mut(coefficients: &mut [BigInt], modulus: &BigInt) {
    let half_modulus = modulus / BigInt::from(2);
    coefficients
        .iter_mut()
        .for_each(|x| *x = reduce_and_center(x, modulus, &half_modulus));
}
fn reduce_and_center_coefficients(coefficients: &mut [BigInt], modulus: &BigInt) -> Vec<BigInt> {
    let half_modulus = modulus / BigInt::from(2);
    coefficients
        .iter()
        .map(|x| reduce_and_center(x, modulus, &half_modulus))
        .collect()
}

/// Reduces a polynomial's coefficients within a polynomial ring defined by a cyclotomic polynomial and a modulus.
///
/// This function performs two reductions on the polynomial represented by `coefficients`:
/// 1. **Cyclotomic Reduction**: Reduces the polynomial by the cyclotomic polynomial, replacing
///    the original coefficients with the remainder after polynomial division.
/// 2. **Modulus Reduction**: Reduces the coefficients of the polynomial modulo a given modulus,
///    centering the coefficients within the range [-modulus/2, modulus/2).
///
/// # Arguments
///
/// * `coefficients` - A mutable reference to a `Vec<BigInt>` representing the coefficients of the polynomial
///   to be reduced. The coefficients should be in descending order of degree.
/// * `cyclo` - A slice of `BigInt` representing the coefficients of the cyclotomic polynomial (typically x^N + 1).
/// * `modulus` - A reference to a `BigInt` representing the modulus for the coefficient reduction. The coefficients
///   will be reduced and centered modulo this value.
fn reduce_in_ring(coefficients: &mut Vec<BigInt>, cyclo: &[BigInt], modulus: &BigInt) {
    reduce_coefficients_by_cyclo(coefficients, cyclo);
    reduce_and_center_coefficients_mut(coefficients, modulus);
}

/// Reduces each element in the given slice of `BigInt` by the modulus `p`.
///
/// This function takes a slice of `BigInt` coefficients and applies the modulus operation
/// on each element. It ensures the result is within the range `[0, p-1]` by adding `p`
/// before applying the modulus operation. The result is collected into a new `Vec<BigInt>`.
///
/// # Arguments
///
/// * `coefficients` - A slice of `BigInt` representing the coefficients to be reduced.
/// * `p` - A reference to a `BigInt` that represents the modulus value.
///
/// # Returns
///
/// A `Vec<BigInt>` where each element is reduced modulo `p`.
fn reduce_coefficients(coefficients: &[BigInt], p: &BigInt) -> Vec<BigInt> {
    coefficients.iter().map(|coeff| (coeff + p) % p).collect()
}

fn reduce_coefficients_2d(coefficient_matrix: &[Vec<BigInt>], p: &BigInt) -> Vec<Vec<BigInt>> {
    coefficient_matrix
        .iter()
        .map(|coeffs| reduce_coefficients(coeffs, p))
        .collect()
}

/// Mutably reduces each element in the given slice of `BigInt` by the modulus `p`.
///
/// This function modifies the given mutable slice of `BigInt` coefficients in place. It adds `p`
/// to each element before applying the modulus operation, ensuring the results are within the range `[0, p-1]`.
///
/// # Arguments
///
/// * `coefficients` - A mutable slice of `BigInt` representing the coefficients to be reduced.
/// * `p` - A reference to a `BigInt` that represents the modulus value.
fn reduce_coefficients_mut(coefficients: &mut [BigInt], p: &BigInt) {
    for coeff in coefficients.iter_mut() {
        *coeff += p;
        *coeff %= p;
    }
}

fn range_check_centered(vec: &[BigInt], lower_bound: &BigInt, upper_bound: &BigInt) -> bool {
    vec.iter()
        .all(|coeff| coeff >= lower_bound && coeff <= upper_bound)
}

fn range_check_standard(vec: &[BigInt], bound: &BigInt, modulus: &BigInt) -> bool {
    vec.iter().all(|coeff| {
        (coeff >= &BigInt::from(0) && coeff <= bound)
            || (coeff >= &(modulus - bound) && coeff < modulus)
    })
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
