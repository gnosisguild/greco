use fhe::bfv::{BfvParametersBuilder, Encoding, Plaintext, PublicKey, SecretKey};
use fhe_math::rq::{Poly, Representation};
use fhe_math::zq::Modulus;
use fhe_traits::FheEncoder;
use itertools::izip;
use ndarray::{s, Array1, Array2, ArrayView2};
use num_bigint::BigInt;
use num_traits::{Num, ToPrimitive, Zero};
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::json;
use std::ops::Deref;
use std::{fs, vec};

fn main() {
    // Set up the BFV parameters
    let N: u64 = 128;
    let plaintext_modulus: u64 = 65537;
    let moduli: Vec<u64> = vec![0xffffffff00001, 0xfffffffe40001];

    let params = BfvParametersBuilder::new()
        .set_degree(N as usize)
        .set_plaintext_modulus(plaintext_modulus)
        .set_moduli(&moduli)
        .build_arc()
        .unwrap();

    // Use a seedable rng for experimental reproducibility
    let mut rng = StdRng::seed_from_u64(0);

    // Generate the secret and public keys
    let sk = SecretKey::random(&params, &mut rng);
    let pk = PublicKey::new(&sk, &mut rng);

    // Encrypt the input plaintext
    let input: Vec<u64> = vec![0, 1];
    let pt = Plaintext::try_encode(&input, Encoding::poly(), &params).unwrap();
    let (ct, mut u_rns, mut e0_rns, mut e1_rns) = pk.try_encrypt_extended(&pt, &mut rng).unwrap();

    // Extract parameters (context and plaintext modulus)
    let ctx = params.ctx_at_level(pt.level()).unwrap().clone();
    let t = Modulus::new(params.plaintext()).unwrap();

    // Calculate k1 (independent of qi), and reverse
    let q_mod_t = (ctx.modulus() % t.modulus()).to_u64().unwrap(); // [q]_t
    let mut k1_u64 = pt.value.deref().to_vec(); // m
    t.scalar_mul_vec(&mut k1_u64, q_mod_t); // k1 = [q*m]_t
    let mut k1: Vec<BigInt> = k1_u64.iter().map(|&x| BigInt::from(x)).rev().collect();

    // Extract single vectors of u, e1, and e2, as Vec<BigInt>, and reverse
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

    // Perform the main computation logic
    let mut r2is = Vec::new();
    let mut r1is = Vec::new();
    let mut p2is = Vec::new();
    let mut p1is = Vec::new();
    for (i, (qi, ct0_coeffs, ct1_coeffs, pk0_coeffs, pk1_coeffs)) in izip!(
        ctx.moduli_operators(),
        ct0.coefficients().outer_iter(),
        ct1.coefficients().outer_iter(),
        pk0.coefficients().outer_iter(),
        pk1.coefficients().outer_iter(),
    )
    .enumerate()
    {
        // Convert to vectors of bigint and reverse order.
        let ct0i: Vec<BigInt> = ct0_coeffs.iter().map(|&x| BigInt::from(x)).rev().collect();
        let ct1i: Vec<BigInt> = ct1_coeffs.iter().map(|&x| BigInt::from(x)).rev().collect();
        let pk0i: Vec<BigInt> = pk0_coeffs.iter().map(|&x| BigInt::from(x)).rev().collect();
        let pk1i: Vec<BigInt> = pk1_coeffs.iter().map(|&x| BigInt::from(x)).rev().collect();

        // k0qi = -t^{-1} mod qi
        let koqi_u64 = qi.inv(qi.neg(t.modulus())).unwrap();
        let k0qi = BigInt::from(koqi_u64);

        // k1 * k0qi as big int
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

        // Compute r2i numerator = ct0i - ct0i_hat
        let ct0i_minus_ct0i_hat = poly_sub(&ct0i, &ct0i_hat);
        assert_eq!((ct0i_minus_ct0i_hat.len() as u64) - 1, 2 * (N - 1));

        // Compute r2i as the quotient of numerator divided by the
        // cyclotomic polynomial, and reduce/center the resulting
        // coefficients to produce: (ct0i - ct0i_hat) / (x^N + 1) mod Z_qi
        let (mut r2i, r2i_rem) = poly_div(&ct0i_minus_ct0i_hat, &cyclo);
        assert_eq!((r2i.len() as u64) - 1, N - 2); // Order(r2i) = N - 2
        assert_eq!(
            ct0i_minus_ct0i_hat,
            poly_add(&poly_mul(&r2i, &cyclo), &r2i_rem)
        );
        reduce_and_center_coefficients(&mut r2i, &BigInt::from(qi.modulus()));

        // Calculate r1i using r2i and cyclo
        let r1i_num: Vec<BigInt> = {
            let r2i_times_cyclo = poly_mul(&r2i, &cyclo);
            assert_eq!((r2i_times_cyclo.len() as u64) - 1, 2 * (N - 1));

            poly_sub(&ct0i_minus_ct0i_hat, &r2i_times_cyclo)
        };
        assert_eq!((r1i_num.len() as u64) - 1, 2 * (N - 1));
        let (r1i, r1i_rem) = poly_div(&r1i_num, &[BigInt::from(qi.modulus())]);
        assert_eq!((r1i.len() as u64) - 1, 2 * (N - 1)); // Order(r1i) = 2*(N-1)
        assert_eq!(
            r1i_num,
            poly_add(&poly_mul(&r1i, &[BigInt::from(qi.modulus())]), &r1i_rem)
        );

        // Calculate ct1i_hat = pk1i * ui + e1i
        let ct1i_hat = {
            let pk1i_times_u = poly_mul(&pk1i, &u);
            assert_eq!((pk1i_times_u.len() as u64) - 1, 2 * (N - 1));

            poly_add(&pk1i_times_u, &e1)
        };
        assert_eq!((ct1i_hat.len() as u64) - 1, 2 * (N - 1));

        // Compute p2i numerator = ct1i - ct1i_hat
        let ct1i_minus_ct1i_hat = poly_sub(&ct1i, &ct1i_hat);
        assert_eq!((ct1i_minus_ct1i_hat.len() as u64) - 1, 2 * (N - 1));

        // Compute p2i as the quotient of numerator divided by the
        // cyclotomic polynomial, and reduce/center the resulting
        // coefficients to produce: (ct1i - ct1i_hat) / (x^N + 1) mod Z_qi
        let (mut p2i, p2i_rem) = poly_div(&ct1i_minus_ct1i_hat, &cyclo.clone());
        assert_eq!((p2i.len() as u64) - 1, N - 2); // Order(p2i) = N - 2
        assert_eq!(
            ct1i_minus_ct1i_hat,
            poly_add(&poly_mul(&p2i, &cyclo), &p2i_rem)
        );
        reduce_and_center_coefficients(&mut p2i, &BigInt::from(moduli[i]));

        // Calculate p1i using p2i and cyclo
        let p1i_num = {
            let p2i_times_cyclo = poly_mul(&p2i, &cyclo);
            assert_eq!((p2i_times_cyclo.len() as u64) - 1, 2 * (N - 1));

            poly_sub(&ct1i_minus_ct1i_hat, &p2i_times_cyclo)
        };
        assert_eq!((p1i_num.len() as u64) - 1, 2 * (N - 1));
        let (p1i, p1i_rem) = poly_div(&p1i_num, &[BigInt::from(qi.modulus())]);
        assert_eq!((p1i.len() as u64) - 1, 2 * (N - 1)); // Order(p1i) = 2*(N-1)
        assert_eq!(
            p1i_num,
            poly_add(&poly_mul(&p1i, &[BigInt::from(qi.modulus())]), &p1i_rem)
        );

        r2is.push(r2i);
        r1is.push(r1i);
        p2is.push(p2i);
        p1is.push(p1i);
    }

    for r2i in r2is.iter_mut() {
        add_p_and_modulo(r2i, &p);
    }
    for r1i in r1is.iter_mut() {
        add_p_and_modulo(r1i, &p);
    }
    for p2i in p2is.iter_mut() {
        add_p_and_modulo(p2i, &p);
    }
    for p1i in p1is.iter_mut() {
        add_p_and_modulo(p1i, &p);
    }
    add_p_and_modulo(&mut k1, &p);

    let pk_array_assigned = add_p_and_modulo_poly(&pk0, &p);
    let pk1_array_assigned = add_p_and_modulo_poly(&pk1, &p);
    let u_assigned = add_p_and_modulo_poly(&u_rns, &p);
    let e0_assigned = add_p_and_modulo_poly(&e0_rns, &p);
    let e1_assigned = add_p_and_modulo_poly(&e1_rns, &p);
    let ct0_assigned = add_p_and_modulo_poly(&ct0, &p);
    let ct1_assigned = add_p_and_modulo_poly(&ct1, &p);

    let json_data = json!({
        "pk0_qi": to_string_2d_vec(&pk_array_assigned),
        "pk1_qi": to_string_2d_vec(&pk1_array_assigned),
        "u": to_string_2d_vec(&u_assigned),
        "e0": to_string_2d_vec(&e0_assigned),
        "e1": to_string_2d_vec(&e1_assigned),
        "k1": to_string_1d_vec(&k1),
        "r2is": to_string_2d_vec(&r2is),
        "r1is": to_string_2d_vec(&r1is),
        "p2is": to_string_2d_vec(&p2is),
        "p1is": to_string_2d_vec(&p1is),
        "ct0is": to_string_2d_vec(&ct0_assigned),
        "ct1is": to_string_2d_vec(&ct1_assigned)
    });

    fs::write("output.json", json_data.to_string()).unwrap();
    fs::write(
        "output_zeroes.json",
        create_zeroes_json(params.degree(), params.moduli().len()).to_string(),
    )
    .unwrap();
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

    while remainder.len() > 1 && remainder[0].is_zero() {
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

/// Reduces and centers polynomial coefficients modulo a prime modulus.
///
/// This function takes a mutable slice of polynomial coefficients and reduces each coefficient
/// modulo a prime modulus, adjusting them to be within the range [−(modulus−1)/2, (modulus−1)/2].
///
/// # Parameters
///
/// - `coefficients`: Mutable slice of `BigInt` coefficients to be centered.
/// - `modulus`: Prime modulus used for reduction and centering.
///
/// # Panics
///
/// - Panics if `modulus` is zero due to division by zero.
fn reduce_and_center_coefficients(coefficients: &mut [BigInt], modulus: &BigInt) {
    let half_modulus = modulus / BigInt::from(2);

    for c in coefficients.iter_mut() {
        // Calculate the remainder
        let mut r = &*c % modulus;

        // Adjust the remainder if it is greater than half_modulus
        if r > half_modulus {
            r -= modulus;
        }

        // Assign the centered coefficient back
        *c = r;
    }
}

fn add_p_and_modulo(coefficients: &mut [BigInt], p: &BigInt) {
    for coeff in coefficients.iter_mut() {
        *coeff += p;
        *coeff %= p;
    }
}

fn add_p_and_modulo_poly(poly: &Poly, p: &BigInt) -> Vec<Vec<BigInt>> {
    poly.coefficients()
        .outer_iter()
        .map(|row| {
            row.iter()
                .map(|&coef| {
                    let mut big_int_coef = BigInt::from(coef);
                    big_int_coef += p;
                    big_int_coef %= p;
                    big_int_coef
                })
                .collect::<Vec<BigInt>>()
        })
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
