use fhe::bfv::{BfvParametersBuilder, Ciphertext, Encoding, Plaintext, PublicKey, SecretKey};
use fhe_math::rns::ScalingFactor;
use fhe_math::rq::{scaler::Scaler, traits::TryConvertFrom, Context, Poly, Representation};
use fhe_traits::{FheEncoder, FheEncrypter};
use num_bigint::{BigInt, BigUint};
use num_traits::{Num, ToPrimitive};
use rand::thread_rng;
use serde_json::json;
use std::ops::Neg;
use std::{fs, vec};
use ndarray::Array2;


// Still WIP
fn poly_div(num: &Poly, cyclo: &Poly) -> (Poly, Poly) {
    let ctx = num.ctx().clone();
    let degree = ctx.degree;

    let mut quotient_coeffs = Array2::zeros((num.coefficients().nrows(), degree));
    let mut remainder_coeffs = num.coefficients().to_owned();

    for i in 0..num.coefficients().nrows() {
        let divisor_leading_coeff = cyclo.coefficients()[(i, 0)];

        for j in 0..(degree - cyclo.coefficients().ncols() + 1) {
            let coeff = remainder_coeffs[(i, j)] / divisor_leading_coeff;
            quotient_coeffs[(i, j)] = coeff;

            for k in 0..cyclo.coefficients().ncols() {
                remainder_coeffs[(i, j + k)] = remainder_coeffs[(i, j + k)]
                    .wrapping_sub(cyclo.coefficients()[(i, k)].wrapping_mul(coeff));
            }
        }
    }

    // Create Poly objects for quotient and remainder
    let quotient_poly = Poly::try_convert_from(
        quotient_coeffs.as_slice().unwrap(),
        &ctx,
        false,
        Representation::PowerBasis,
    )
    .unwrap();

    let remainder_poly = Poly::try_convert_from(
        remainder_coeffs.as_slice().unwrap(),
        &ctx,
        false,
        Representation::PowerBasis,
    )
    .unwrap();

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
    // TODO: return u, e0, e1 from pk.try_encrypt
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
        &params.ctx_at_level(pt.level()).unwrap(),
        &params.ctx_at_level(pt.level()).unwrap(),
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

    let mut r2is = Vec::new();
    let mut r1is = Vec::new();
    let mut k0is = Vec::new();
    let mut ct0is = Vec::new();
    let mut ct0is_hat = Vec::new();
    let mut ct1is = Vec::new();
    let mut ct1is_hat = Vec::new();
    let mut pk0is = Vec::new();
    let mut pk1is = Vec::new();
    let mut p1is = Vec::new();
    let mut p2is = Vec::new();

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

    let ct0 = &ct.c[0]; // First polynomial (corresponding to ct0)
    let ct1 = &ct.c[1]; // Second polynomial (corresponding to ct1)
    let pk0 = &pk.c.c[0]; // First polynomial (corresponding to pk0)
    let pk1 = &pk.c.c[1]; // Second polynomial (corresponding to pk1)

    for (modulus_index, (ct0_coeffs, ct1_coeffs)) in ct0
        .coefficients()
        .outer_iter()
        .zip(ct1.coefficients().outer_iter())
        .enumerate()
    {
        // ct0i = cti[0]
        let mut ct0i = Poly::try_convert_from(
            &ct0_coeffs.iter().cloned().collect::<Vec<u64>>(),
            &params.ctx_at_level(modulus_index).unwrap(),
            false,
            Representation::PowerBasis,
        )
        .unwrap();
        ct0i.change_representation(Representation::Ntt);

        // ct1i = cti[1]
        let mut ct1i = Poly::try_convert_from(
            &ct1_coeffs.iter().cloned().collect::<Vec<u64>>(),
            &params.ctx_at_level(modulus_index).unwrap(),
            false,
            Representation::PowerBasis,
        )
        .unwrap();
        ct1i.change_representation(Representation::Ntt);

        // pk0 = Polynomial(pk_array[i])
        let mut pk0i = Poly::try_convert_from(
            &pk0.coefficients()
                .outer_iter()
                .nth(modulus_index)
                .unwrap()
                .iter()
                .cloned()
                .collect::<Vec<u64>>(),
            &params.ctx_at_level(modulus_index).unwrap(),
            false,
            Representation::PowerBasis,
        )
        .unwrap();
        pk0i.change_representation(Representation::Ntt);

        // pk1 = Polynomial(pk_array[i])
        let mut pk1i = Poly::try_convert_from(
            &pk1.coefficients()
                .outer_iter()
                .nth(modulus_index)
                .unwrap()
                .iter()
                .cloned()
                .collect::<Vec<u64>>(),
            &params.ctx_at_level(modulus_index).unwrap(),
            false,
            Representation::PowerBasis,
        )
        .unwrap();
        pk1i.change_representation(Representation::Ntt);

        // k0i = pow(-t, -1, qis[i])
        let k0i = BigInt::from(plaintext_modulus)
            .neg()
            .modinv(&BigInt::from(moduli[modulus_index]))
            .unwrap()
            .to_u64()
            .unwrap();

        let scaling_factor = ScalingFactor::new(&BigUint::from(k0i), &BigUint::from(1u64));
        let scaler = Scaler::new(
            &params.ctx_at_level(modulus_index).unwrap(),
            &params.ctx_at_level(modulus_index).unwrap(),
            scaling_factor,
        )
        .unwrap();

        let ct0i_hat: Poly = &pk0i * &u + e0.clone() + k1.scale(&scaler).unwrap();
        let mut num: Poly = ct0i.clone() + ct0i_hat.clone().neg();
        let (mut r2i, _remainder) = poly_div(&num, &cyclo);
        r2i.change_representation(Representation::Ntt);

        num = num + (&r2i * &cyclo).neg();
        let (mut r1i, _remainder) = poly_div(&num, &cyclo);
        r1i.change_representation(Representation::Ntt);

        let ct1i_hat: Poly = &pk1i * &u + e1.clone();
        let num: Poly = ct1i.clone() + ct1i_hat.clone().neg();
        let (mut p2i, _remainder) = poly_div(&num, &cyclo);
        p2i.change_representation(Representation::Ntt);

        let num: Poly = num + (&p2i * &cyclo).neg();
        let (mut p1i, _remainder) = poly_div(&num, &cyclo);
        p1i.change_representation(Representation::Ntt);

        pk1is.push(pk1i);
        p2is.push(p2i);
        p1is.push(p1i);
        ct1is.push(ct1i);
        ct1is_hat.push(ct1i_hat);
        r2is.push(r2i);
        r1is.push(r1i);
        k0is.push(k0i);
        ct0is.push(ct0i);
        ct0is_hat.push(ct0i_hat);
        pk0is.push(pk0i);
    }

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
