use fhe::bfv::{BfvParametersBuilder, Ciphertext, Encoding, Plaintext, PublicKey, SecretKey};
use fhe_math::rns::ScalingFactor;
use fhe_math::rq::{scaler::Scaler, traits::TryConvertFrom, Context, Poly, Representation};
use fhe_traits::{FheEncoder, FheEncrypter};
use ndarray::Array2;
use num_bigint::{BigInt, BigUint};
use num_traits::{Num, ToPrimitive};
use rand::thread_rng;
use serde_json::json;
use std::ops::Neg;
use std::{fs, vec};

fn poly_div(num: &Poly, divisor: &Poly, modulus: Option<u64>) -> (Poly, Poly) {
    let ctx = num.ctx().clone();
    let degree = ctx.degree;

    let mut quotient_coeffs = Array2::zeros((num.coefficients().nrows(), degree));
    let mut remainder_coeffs = num.coefficients().to_owned();

    let is_div_by_modulus = modulus.is_some();

    for i in 0..num.coefficients().nrows() {
        if is_div_by_modulus {
            let modulus_value = modulus.unwrap();

            for j in 0..degree {
                quotient_coeffs[(i, j)] = remainder_coeffs[(i, j)] / modulus_value;
                remainder_coeffs[(i, j)] %= modulus_value;
            }
        } else {
            let divisor_leading_coeff = divisor.coefficients()[(i, 0)];

            for j in 0..(degree - divisor.coefficients().ncols() + 1) {
                let coeff = remainder_coeffs[(i, j)] / divisor_leading_coeff;
                quotient_coeffs[(i, j)] = coeff;

                for k in 0..divisor.coefficients().ncols() {
                    remainder_coeffs[(i, j + k)] = remainder_coeffs[(i, j + k)]
                        .wrapping_sub(divisor.coefficients()[(i, k)].wrapping_mul(coeff));
                }
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
    let degree: u64 = 128;
    let plaintext_modulus: u64 = 65537;
    let moduli: Vec<u64> = vec![0xffffffff00001, 0xfffffffe40001];

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
    let (ct, u, e0, e1) = pk.try_encrypt_extended(&pt, &mut rng).unwrap();

    // Context for the plaintext 
    let ctx = params.ctx_at_level(pt.level()).unwrap();


    // k1 = m.scalar_mul(crt_moduli.q)
    let q: BigUint = ctx.rns.product.clone();
    let scaling_factor = ScalingFactor::new(&q, &BigUint::from(1u64));
    let scaler = Scaler::new(ctx, ctx, scaling_factor).unwrap();
    let mut k1 = pt.to_poly().scale(&scaler).unwrap();
    k1.change_representation(Representation::Ntt);

    let p = BigUint::from_str_radix(
        "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        10,
    )
    .unwrap();

    let mut r2is = Vec::new();
    let mut r1is = Vec::new();
    let mut ct0is = Vec::new();
    let mut ct1is = Vec::new();
    let mut pk0is = Vec::new();
    let mut pk1is = Vec::new();
    let mut p1is = Vec::new();
    let mut p2is = Vec::new();

    let mut ct0: Poly = ct.c[0].clone();
    ct0.change_representation(Representation::PowerBasis);
    let mut ct1: Poly = ct.c[1].clone();
    ct1.change_representation(Representation::PowerBasis);
    let mut pk_array: Poly = pk.c.c[0].clone();
    pk_array.change_representation(Representation::PowerBasis);
    let mut pk1_array: Poly = pk.c.c[1].clone();
    pk1_array.change_representation(Representation::PowerBasis);

    let mut cyclo_coeffs = vec![0u64; degree as usize];
    cyclo_coeffs[0] = 1; // x^0 term
    cyclo_coeffs[(degree - 1) as usize] = 1; // x^(n-1) term
    let mut cyclo =
        Poly::try_convert_from(cyclo_coeffs, ctx, false, Representation::PowerBasis).unwrap();

    for (modulus_index, (ct0_coeffs, ct1_coeffs)) in ct0
        .coefficients()
        .outer_iter()
        .zip(ct1.coefficients().outer_iter())
        .enumerate()
    {
        // ct0i = cti[0]
        let mut ct0i = Poly::try_convert_from(
            ct0_coeffs.as_slice().unwrap(),
            ctx,
            false,
            Representation::PowerBasis,
        )
        .unwrap();

        // ct1i = cti[1]
        let mut ct1i = Poly::try_convert_from(
            ct1_coeffs.as_slice().unwrap(),
            ctx,
            false,
            Representation::PowerBasis,
        )
        .unwrap();

        // pk0 = Polynomial(pk_array[i])
        let mut pk0 = Poly::try_convert_from(
            pk_array
                .coefficients()
                .outer_iter()
                .nth(modulus_index)
                .unwrap()
                .as_slice()
                .unwrap(),
            ctx,
            false,
            Representation::PowerBasis,
        )
        .unwrap();

        // pk1 = Polynomial(pk_array[i])
        let mut pk1 = Poly::try_convert_from(
            pk1_array
                .coefficients()
                .outer_iter()
                .nth(modulus_index)
                .unwrap()
                .as_slice()
                .unwrap(),
            ctx,
            false,
            Representation::PowerBasis,
        )
        .unwrap();

        // k0i = pow(-t, -1, qis[i])
        let k0i: u64 = BigInt::from(plaintext_modulus)
            .neg()
            .modinv(&BigInt::from(moduli[modulus_index]))
            .unwrap()
            .to_u64()
            .unwrap();

        let scaling_factor = ScalingFactor::new(&BigUint::from(k0i), &BigUint::from(1u64));
        let scaler = Scaler::new(ctx, ctx, scaling_factor).unwrap();

        pk0.change_representation(Representation::Ntt);
        pk1.change_representation(Representation::Ntt);

        // ct0i_hat = pk0 * u + e0 + k1.scalar_mul(k0i)
        let mut ct0i_hat: Poly = &pk0 * &u + e0.clone() + k1.scale(&scaler).unwrap();
        ct0i_hat.change_representation(Representation::PowerBasis);
        // num = ct0i + ct0i_hat.scalar_mul(-1)
        let mut num: Poly = ct0i.clone() + ct0i_hat.clone().neg();

        let (mut r2i, _remainder) = poly_div(&num, &cyclo, None);
        r2i.change_representation(Representation::Ntt);
        cyclo.change_representation(Representation::Ntt);
        num.change_representation(Representation::Ntt);

        num = num + (&r2i * &cyclo).neg();
        num.change_representation(Representation::PowerBasis);
        let (mut r1i, _remainder) = poly_div(&num, &cyclo, Some(moduli[modulus_index]));
        r1i.change_representation(Representation::Ntt);

        // ct1i_hat = pk1 * u + e1
        let mut ct1i_hat: Poly = &pk1 * &u + e1.clone();
        ct1i_hat.change_representation(Representation::PowerBasis);
        // num = ct1i + ct1i_hat.scalar_mul(-1)
        let mut num: Poly = ct1i.clone() + ct1i_hat.clone().neg();
        num.change_representation(Representation::PowerBasis);
        cyclo.change_representation(Representation::PowerBasis);
        let (mut p2i, _remainder) = poly_div(&num, &cyclo, None);
        p2i.change_representation(Representation::Ntt);

        num.change_representation(Representation::Ntt);
        cyclo.change_representation(Representation::Ntt);
        let mut num: Poly = num + (&p2i * &cyclo).neg();
        num.change_representation(Representation::PowerBasis);
        let (mut p1i, _remainder) = poly_div(&num, &&cyclo, Some(moduli[modulus_index]));
        p1i.change_representation(Representation::Ntt);

        pk1.change_representation(Representation::PowerBasis);
        pk1is.push(pk1);
        p2i.change_representation(Representation::PowerBasis);
        p2is.push(p2i);
        p1i.change_representation(Representation::PowerBasis);
        p1is.push(p1i);
        ct1i.change_representation(Representation::PowerBasis);
        ct1is.push(ct1i);
        r2i.change_representation(Representation::PowerBasis);
        r2is.push(r2i);
        r1i.change_representation(Representation::PowerBasis);
        r1is.push(r1i);
        ct0i.change_representation(Representation::PowerBasis);
        ct0is.push(ct0i);
        pk0.change_representation(Representation::PowerBasis);
        pk0is.push(pk0);
    }


    let pk0_qi: Vec<Vec<String>> = pk_array
        .coefficients()
        .outer_iter()
        .map(|row| row.iter().map(|coef| coef.to_string()).collect())
        .collect();
    let pk1_qi: Vec<Vec<String>> = pk1_array
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

    let ct0is_str: Vec<Vec<String>> = ct0
        .coefficients()
        .outer_iter()
        .map(|row| row.iter().map(|coef| coef.to_string()).collect())
        .collect();

    let ct1is_str: Vec<Vec<String>> = ct1
        .coefficients()
        .outer_iter()
        .map(|row| row.iter().map(|coef| coef.to_string()).collect())
        .collect();

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
        "ct1is": ct1is_str
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
        "ct1is": vec![vec![String::from("0"); degree as usize]; moduli.len()]
    });

    let json_zeroes_output_path = "output_zeroes.json";
    fs::write(json_zeroes_output_path, json_data_zeroes.to_string()).unwrap();
}
