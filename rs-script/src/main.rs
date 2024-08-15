use fhe::bfv::{BfvParametersBuilder, Ciphertext, Encoding, Plaintext, PublicKey, SecretKey};
use fhe_math::rns::ScalingFactor;
use fhe_math::rq::{scaler::Scaler, traits::TryConvertFrom, Context, Poly, Representation};
use fhe_traits::{FheEncoder, FheEncrypter};
use num_bigint::BigUint;
use num_traits::{Num, ToPrimitive};
use rand::thread_rng;
use serde_json::json;
use std::{fs, vec};
use std::ops::Neg;
use std::sync::Arc;

fn poly_div(num: &Poly, cyclo: &Poly, ctx: &Arc<Context>) -> (Poly, Poly) {
    let num_coeffs: Vec<u64> = num.coefficients().iter().cloned().collect();
    let cyclo_coeffs: Vec<u64> = cyclo.coefficients().iter().cloned().collect();

    let mut quotient = vec![0u64; num_coeffs.len() - cyclo_coeffs.len() + 1];
    let mut remainder = num_coeffs.clone();

    let divisor_leading_coeff = cyclo_coeffs[0];

    for i in 0..quotient.len() {
        let coeff = remainder[i] / divisor_leading_coeff;
        quotient[i] = coeff;

        for j in 0..cyclo_coeffs.len() {
            remainder[i + j] = remainder[i + j].wrapping_sub(cyclo_coeffs[j].wrapping_mul(coeff));
        }
    }

    let quotient_poly =
        Poly::try_convert_from(&quotient, ctx, false, Representation::PowerBasis).unwrap();
    let remainder_poly =
        Poly::try_convert_from(&remainder, ctx, false, Representation::PowerBasis).unwrap();

    (quotient_poly, remainder_poly)
}

fn main() {
    let degree: u64 = 1024;
    let plaintext_modulus: u64 = 65537;
    let moduli: Vec<u64> = vec![
        1152921504606584833,
        1152921504598720513,
        1152921504597016577,
        1152921504595968001,
        1152921504595640321,
        1152921504593412097,
        1152921504592822273,
        1152921504592429057,
        1152921504589938689,
        1152921504586530817,
        1152921504585547777,
        1152921504583647233,
        1152921504581877761,
        1152921504581419009,
        1152921504580894721,
    ];

    let params = BfvParametersBuilder::new()
        .set_degree(degree as usize)
        .set_plaintext_modulus(plaintext_modulus)
        .set_moduli(&moduli)
        .build_arc()
        .unwrap();

    let mut rng = thread_rng();

    let sk = SecretKey::random(&params, &mut rng);
    let pk = PublicKey::new(&sk, &mut rng);

    let input: Vec<u64> = vec![0, 1];
    let pt = Plaintext::try_encode(&input, Encoding::poly(), &params).unwrap();
    let ct: Ciphertext = pk.try_encrypt(&pt, &mut rng).unwrap();

    let ctx = params.ctx().get(0).unwrap();
    // TODO: Check if u, e0, e1 are the same as used in the pk encryption
    let u = Poly::small(ctx, Representation::Ntt, params.variance(), &mut rng).unwrap();
    let e0 = Poly::small(ctx, Representation::Ntt, params.variance(), &mut rng).unwrap();
    let e1 = Poly::small(ctx, Representation::Ntt, params.variance(), &mut rng).unwrap();

    // Scaling the message polynomial to compute k1
    let q: BigUint = params
        .moduli()
        .iter()
        .map(|modulus| BigUint::from(*modulus))
        .product();
    let scaling_factor = ScalingFactor::new(&q, &BigUint::from(1u64));
    let scaler = Scaler::new(
        &params.ctx_at_level(0).unwrap(),
        &params.ctx_at_level(0).unwrap(),
        scaling_factor,
    )
    .unwrap();

    let mut k1 = pt.to_poly().scale(&scaler).unwrap();
    k1.change_representation(Representation::Ntt);

    let p = BigUint::from_str_radix(
        "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        10,
    )
    .unwrap();

    let pk0_array: Vec<u64> = pk.c.c[0].coefficients().iter().copied().collect();
    let pk1_array: Vec<u64> = pk.c.c[1].coefficients().iter().copied().collect();

    let mut cyclo_coeffs = vec![0u64; degree as usize];
    cyclo_coeffs[0] = 1; // x^0 term
    cyclo_coeffs[(degree - 1) as usize] = 1; // x^(n-1) term
    let mut cyclo = Poly::try_convert_from(
        &cyclo_coeffs,
        &params.ctx_at_level(0).unwrap(),
        false,
        Representation::PowerBasis,
    )
    .unwrap();

    cyclo.change_representation(Representation::Ntt);

    let ct0i = &ct.c[0];
    let ct1i = &ct.c[1];

    let pk0 = Poly::try_convert_from(&pk0_array, &ctx, false, Representation::Ntt).unwrap();
    let pk1 = Poly::try_convert_from(&pk1_array, &ctx, false, Representation::Ntt).unwrap();

    let k0i = BigUint::from(plaintext_modulus)
    .modinv(&BigUint::from(moduli[0]))
        .unwrap()
        .to_u64()
        .unwrap();

    let scaling_factor = ScalingFactor::new(&BigUint::from(k0i), &BigUint::from(1u64));
    let scaler = Scaler::new(
        &params.ctx_at_level(0).unwrap(),
        &params.ctx_at_level(0).unwrap(),
        scaling_factor,
    )
    .unwrap();

    println!("pk0: {:?}", pk0.representation());
    println!("pk1: {:?}", pk1.representation());
    println!("u: {:?}", u.representation());
    println!("e0: {:?}", e0.representation());
    println!("e1: {:?}", e1.representation());
    println!("cyclo: {:?}", cyclo.representation());

    let ct0_hat = &pk0 * &u + e0.clone() + k1.scale(&scaler).unwrap();
    let mut num: Poly = ct0i - &ct0_hat;
    let (mut r2, _remainder) = poly_div(&num, &cyclo, &ctx);
    r2.change_representation(Representation::Ntt);

    num = num + (&r2 * &cyclo).neg();

    let q_coeffs = q.to_u64_digits();
    let q_poly = Poly::try_convert_from(
        &q_coeffs,
        &params.ctx_at_level(0).unwrap(),
        false,
        Representation::PowerBasis,
    )
    .unwrap();
    let (mut r1, _remainder) = poly_div(&num, &q_poly, &ctx);
    r1.change_representation(Representation::Ntt);

    let mut ct1_hat = &pk1 * &u;
    ct1_hat = &ct1_hat * &e1;

    let num = ct1i - &ct1_hat;
    let (mut p2, _remainder) = poly_div(&num, &cyclo, &ctx);
    p2.change_representation(Representation::Ntt);

    let num = num + (&p2 * &cyclo).neg();
    let (mut p1, _remainder) = poly_div(&num, &q_poly, &ctx);
    p1.change_representation(Representation::Ntt);

    let r2is = vec![r2];
    let r1is = vec![r1];
    let k0is =  vec![k0i];
    let ct0is = vec![ct0i.clone()];
    let ct0is_hat = vec![ct0_hat];
    let ct1is = vec![ct1i.clone()];
    let ct1is_hat = vec![ct1_hat];
    let pk0is = vec![pk0.clone()];
    let pk1is = vec![pk1.clone()];
    let p1is = vec![p1];
    let p2is = vec![p2];

    let pk0_qi: Vec<Vec<String>> = pk0
        .coefficients()
        .outer_iter()
        .map(|row| row.iter().map(|coef| coef.to_string()).collect())
        .collect();

    let pk1_qi: Vec<Vec<String>> = pk1
        .coefficients()
        .outer_iter()
        .map(|row| row.iter().map(|coef| coef.to_string()).collect())
        .collect();

    let u_str: Vec<String> = u
        .coefficients()
        .outer_iter()
        .map(|row| row.iter().map(|coef| coef.to_string()).collect())
        .collect();

    let e0_str: Vec<String> = e0
        .coefficients()
        .outer_iter()
        .map(|row| row.iter().map(|coef| coef.to_string()).collect())
        .collect();

    let e1_str: Vec<String> = e1
        .coefficients()
        .outer_iter()
        .map(|row| row.iter().map(|coef| coef.to_string()).collect())
        .collect();

    let k1_str: Vec<String> = k1
        .coefficients()
        .outer_iter()
        .map(|row| row.iter().map(|coef| coef.to_string()).collect())
        .collect();

    let r2is_str: Vec<Vec<String>> = r2is
        .iter()
        .map(|r2i: &Poly| {
            r2i.coefficients()
                .outer_iter()
                .map(|row| row.iter().map(|coef| coef.to_string()).collect())
                .collect()
        })
        .collect();

    let r1is_str: Vec<Vec<String>> = r1is
        .iter()
        .map(|r1i: &Poly| {
            r1i.coefficients()
                .outer_iter()
                .map(|row| row.iter().map(|coef| coef.to_string()).collect())
                .collect()
        })
        .collect();

    let p2is_str: Vec<Vec<String>> = p2is
        .iter()
        .map(|p2i: &Poly| {
            p2i.coefficients()
                .outer_iter()
                .map(|row| row.iter().map(|coef| coef.to_string()).collect())
                .collect()
        })
        .collect();

    let p1is_str: Vec<Vec<String>> = p1is
        .iter()
        .map(|p1i: &Poly| {
            p1i.coefficients()
                .outer_iter()
                .map(|row| row.iter().map(|coef| coef.to_string()).collect())
                .collect()
        })
        .collect();

    let ct0is_str: Vec<Vec<String>> = ct0is
        .iter()
        .map(|ct0i: &Poly| {
            ct0i.coefficients()
                .outer_iter()
                .map(|row| row.iter().map(|coef| coef.to_string()).collect())
                .collect()
        })
        .collect();

    let ct1is_str: Vec<Vec<String>> = ct1is
        .iter()
        .map(|ct1i: &Poly| {
            ct1i.coefficients()
                .outer_iter()
                .map(|row| row.iter().map(|coef| coef.to_string()).collect())
                .collect()
        })
        .collect();

    let ct0is_hat_str: Vec<Vec<String>> = ct0is_hat
        .iter()
        .map(|ct0i_hat: &Poly| {
            ct0i_hat
                .coefficients()
                .outer_iter()
                .map(|row| row.iter().map(|coef| coef.to_string()).collect())
                .collect()
        })
        .collect();

    let ct1is_hat_str: Vec<Vec<String>> = ct1is_hat
        .iter()
        .map(|ct1i_hat: &Poly| {
            ct1i_hat
                .coefficients()
                .outer_iter()
                .map(|row| row.iter().map(|coef| coef.to_string()).collect())
                .collect()
        })
        .collect();

    let pk0is_str: Vec<Vec<String>> = pk0is
        .iter()
        .map(|pk0i: &Poly| {
            pk0i.coefficients()
                .outer_iter()
                .map(|row| row.iter().map(|coef| coef.to_string()).collect())
                .collect()
        })
        .collect();

    let pk1is_str: Vec<Vec<String>> = pk1is
        .iter()
        .map(|pk1i: &Poly| {
            pk1i.coefficients()
                .outer_iter()
                .map(|row| row.iter().map(|coef| coef.to_string()).collect())
                .collect()
        })
        .collect();

    let k0is_str: Vec<String> = k0is.iter().map(|coef| coef.to_string()).collect();

    let json_data = json!({
        "pk0_qi": pk0_qi,
        "pk1_qi": pk1_qi,
        "u": u_str,
        "e0": e0_str,
        "e1": e1_str,
        "k1": k1_str,
        "r2is": r2is_str,
        "r1is": r1is_str,
        "p2is": p2is_str,
        "p1is": p1is_str,
        "ct0is": ct0is_str,
        "ct1is": ct1is_str,
        "ct0is_hat": ct0is_hat_str,
        "ct1is_hat": ct1is_hat_str,
        "pk0is": pk0is_str,
        "pk1is": pk1is_str,
        "k0is": k0is_str
    });

    let json_output_path = "output.json";
    fs::write(json_output_path, json_data.to_string()).unwrap();

    let json_data_zeroes = json!({
        "pk0_qi": vec![vec![String::from("0"); degree as usize]; moduli.len()],
        "pk1_qi": vec![vec![String::from("0"); degree as usize]; moduli.len()],
        "u": vec![String::from("0"); degree as usize],
        "e0": vec![String::from("0"); degree as usize],
        "e1": vec![String::from("0"); degree as usize],
        "k1": vec![String::from("0"); degree as usize],
        "r2is": vec![vec![String::from("0"); degree as usize]; moduli.len()],
        "r1is": vec![vec![String::from("0"); degree as usize]; moduli.len()],
        "p2is": vec![vec![String::from("0"); degree as usize]; moduli.len()],
        "p1is": vec![vec![String::from("0"); degree as usize]; moduli.len()],
        "ct0is": vec![vec![String::from("0"); degree as usize]; moduli.len()],
        "ct1is": vec![vec![String::from("0"); degree as usize]; moduli.len()],
        "ct0is_hat": vec![vec![String::from("0"); degree as usize]; moduli.len()],
        "ct1is_hat": vec![vec![String::from("0"); degree as usize]; moduli.len()],
        "pk0is": vec![vec![String::from("0"); degree as usize]; moduli.len()],
        "pk1is": vec![vec![String::from("0"); degree as usize]; moduli.len()],
        "k0is": vec![String::from("0"); moduli.len()]
    });

    let json_zeroes_output_path = "output_zeroes.json";
    fs::write(json_zeroes_output_path, json_data_zeroes.to_string()).unwrap();
}
