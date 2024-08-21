use fhe::bfv::{BfvParametersBuilder, Encoding, Plaintext, PublicKey, SecretKey};
use fhe_math::rq::{ Poly, Representation};
use fhe_traits::FheEncoder;
use num_bigint::{BigInt, BigUint};
use num_traits::{Num, ToPrimitive, Zero};
use rand::thread_rng;
use serde_json::json;
use std::ops::Neg;
use std::{fs, vec};

fn main() {
    // Set up the BFV parameters
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

    // Generate the secret and public keys
    let sk = SecretKey::random(&params, &mut rng);
    let pk = PublicKey::new(&sk, &mut rng);

    // Encrypt the input plaintext
    let input: Vec<u64> = vec![0, 1];
    let pt = Plaintext::try_encode(&input, Encoding::poly(), &params).unwrap();
    let (ct, mut u, mut e0, mut e1) = pk.try_encrypt_extended(&pt, &mut rng).unwrap();

    // Convert the representations of the ciphertexts and key components
    u.change_representation(Representation::PowerBasis);
    e0.change_representation(Representation::PowerBasis);
    e1.change_representation(Representation::PowerBasis);

    // Perform the required computations
    let ctx = params.ctx_at_level(pt.level()).unwrap();
    let q = BigInt::from(ctx.rns.product.clone());
    let pt_bigint: Vec<BigInt> = pt.value.iter().map(|&x| BigInt::from(x)).collect();
    let mut k1: Vec<BigInt> = scalar_mul(&pt_bigint, &q);
    reduce_and_center_coefficients(&mut k1, &BigInt::from(plaintext_modulus));

    let p = BigUint::from_str_radix(
        "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        10,
    )
    .unwrap();

    let mut r2is = Vec::new();
    let mut r1is = Vec::new();
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

    let mut cyclo = vec![BigInt::from(0u64); (degree + 1) as usize];
    cyclo[0] = BigInt::from(1u64); // x^0 term
    cyclo[degree as usize] = BigInt::from(1u64); // x^(n-1) term

    // Perform the main computation logic
    for (modulus_index, (ct0_coeffs, ct1_coeffs)) in ct0
        .coefficients()
        .outer_iter()
        .zip(ct1.coefficients().outer_iter())
        .enumerate()
    {
        // Compute the required values for each modulus
        let ct0i: Vec<BigInt> = ct0_coeffs.iter().map(|&x| BigInt::from(x)).collect();
        let ct1i: Vec<BigInt> = ct1_coeffs.iter().map(|&x| BigInt::from(x)).collect();
        let pk0: Vec<BigInt> = pk_array
            .coefficients()
            .outer_iter()
            .nth(modulus_index)
            .unwrap()
            .iter()
            .map(|&x| BigInt::from(x))
            .collect();
        let pk1: Vec<BigInt> = pk1_array
            .coefficients()
            .outer_iter()
            .nth(modulus_index)
            .unwrap()
            .iter()
            .map(|&x| BigInt::from(x))
            .collect();
        let ui: Vec<BigInt> = u
            .coefficients()
            .outer_iter()
            .nth(modulus_index)
            .unwrap()
            .iter()
            .map(|&x| BigInt::from(x))
            .collect();
        let e0i: Vec<BigInt> = e0
            .coefficients()
            .outer_iter()
            .nth(modulus_index)
            .unwrap()
            .iter()
            .map(|&x| BigInt::from(x))
            .collect();
        let e1i: Vec<BigInt> = e1
            .coefficients()
            .outer_iter()
            .nth(modulus_index)
            .unwrap()
            .iter()
            .map(|&x| BigInt::from(x))
            .collect();

        let k0i: u64 = BigInt::from(plaintext_modulus)
            .neg()
            .modinv(&BigInt::from(moduli[modulus_index]))
            .unwrap()
            .to_u64()
            .unwrap();
        let k1i: Vec<BigInt> = scalar_mul(k1.as_slice(), &BigInt::from(k0i));

        let ct0i_hat = poly_add(poly_mul(pk0, ui.clone()), poly_add(e0i, k1i));
        let mut num = poly_sub(ct0i, ct0i_hat);
        reduce_and_center_coefficients(&mut num, &BigInt::from(moduli[modulus_index]));
        let (r2i, _remainder) = poly_div(&num, &cyclo.clone());

        num = poly_sub(num, poly_mul(r2i.clone(), cyclo.clone()));
        let (r1i, _remainder) = poly_div(&num, &vec![BigInt::from(moduli[modulus_index])]);

        let ct1i_hat = poly_add(poly_mul(pk1, ui), e1i);
        num = poly_sub(ct1i, ct1i_hat);
        reduce_and_center_coefficients(&mut num, &BigInt::from(moduli[modulus_index]));
        let (p2i, _remainder) = poly_div(&num, &cyclo.clone());
        num = poly_sub(num, poly_mul(p2i.clone(), cyclo.clone()));
        let (p1i, _remainder) = poly_div(&num, &vec![BigInt::from(moduli[modulus_index])]);

        p2is.push(p2i);
        p1is.push(p1i);
        r2is.push(r2i);
        r1is.push(r1i);
    }
    let json_data = json!({
        "pk0_qi": to_string_poly(&pk_array),
        "pk1_qi": to_string_poly(&pk1_array),
        "u": to_string_poly(&u),
        "e0": to_string_poly(&e0),
        "e1": to_string_poly(&e1),
        "k1": to_string_1d_vec(&k1),
        "r2is": to_string_2d_vec(&r2is),
        "r1is": to_string_2d_vec(&r1is),
        "p2is": to_string_2d_vec(&p2is),
        "p1is": to_string_2d_vec(&p1is),
        "ct0is": to_string_poly(&ct0),
        "ct1is": to_string_poly(&ct1)
    });

    fs::write("output.json", json_data.to_string()).unwrap();
    fs::write(
        "output_zeroes.json",
        create_zeroes_json(params.degree(), params.moduli().len())
            .to_string()).unwrap();
}

fn poly_add(poly1: Vec<BigInt>, poly2: Vec<BigInt>) -> Vec<BigInt> {
    let max_length = std::cmp::max(poly1.len(), poly2.len());

    let mut result = vec![BigInt::zero(); max_length];

    for i in 0..max_length {
        let p1 = if i < poly1.len() {
            &poly1[i]
        } else {
            &BigInt::zero()
        };
        let p2 = if i < poly2.len() {
            &poly2[i]
        } else {
            &BigInt::zero()
        };
        result[i] = p1 + p2;
    }

    result
}

fn poly_sub(poly1: Vec<BigInt>, poly2: Vec<BigInt>) -> Vec<BigInt> {
    let max_length = std::cmp::max(poly1.len(), poly2.len());

    let mut result = vec![BigInt::zero(); max_length];

    for i in 0..max_length {
        let p1 = if i < poly1.len() {
            &poly1[i]
        } else {
            &BigInt::zero()
        };
        let p2 = if i < poly2.len() {
            &poly2[i]
        } else {
            &BigInt::zero()
        };
        result[i] = p1 - p2;
    }

    result
}

fn poly_mul(poly1: Vec<BigInt>, poly2: Vec<BigInt>) -> Vec<BigInt> {
    let product_len = poly1.len() + poly2.len() - 1;
    let mut product = vec![BigInt::zero(); product_len];

    for i in 0..poly1.len() {
        for j in 0..poly2.len() {
            product[i + j] += &poly1[i] * &poly2[j];
        }
    }

    product
}
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
            let rem = &remainder[i + j] - &divisor[j] * &coeff;
            remainder[i + j] = rem;
        }
    }

    while remainder.len() > 1 && remainder[0].is_zero() {
        remainder.remove(0);
    }

    (quotient, remainder)
}

fn scalar_mul(poly: &[BigInt], scalar: &BigInt) -> Vec<BigInt> {
    poly.iter().map(|coeff| coeff * scalar).collect()
}
fn reduce_and_center_coefficients(coefficients: &mut [BigInt], modulus: &BigInt) {
    let half_modulus = modulus / 2 as u64;

    for coeff in coefficients.iter_mut() {
        let mut r = &*coeff % modulus;
        if r > half_modulus {
            r -= modulus;
        }
        *coeff = r;
    }
}

fn to_string_poly(poly: &Poly) -> Vec<Vec<String>> {
    poly.coefficients()
        .outer_iter()
        .map(|row| row.iter().map(|coef| coef.to_string()).collect())
        .collect()
}
fn to_string_1d_vec(poly: &Vec<BigInt>) -> Vec<String> {
    poly.iter().map(|coef| coef.to_string()).collect()
}

fn to_string_2d_vec(poly: &Vec<Vec<BigInt>>) -> Vec<Vec<String>> {
    poly.iter().map(|row| to_string_1d_vec(row)).collect()
}

fn create_zeroes_json(degree: usize, moduli_len: usize) -> serde_json::Value {
    json!({
        "pk0_qi": vec![vec![String::from("0"); degree]; moduli_len],
        "pk1_qi": vec![vec![String::from("0"); degree]; moduli_len],
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
