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

    // Calculate k1 (independent of qi)
    let q_mod_t = (ctx.modulus() % t.modulus()).to_u64().unwrap(); // [q]_t
    let mut k1_u64 = pt.value.deref().to_vec(); // m
    t.scalar_mul_vec(&mut k1_u64, q_mod_t); // k1 = [q*m]_t
    let mut k1: Vec<BigInt> = k1_u64.iter().map(|&x| BigInt::from(x)).collect(); // for now

    // Extract single vectors of u, e1, and e2, as Vec<BigInt>
    u_rns.change_representation(Representation::PowerBasis);
    e0_rns.change_representation(Representation::PowerBasis);
    e1_rns.change_representation(Representation::PowerBasis);
    let u: Vec<BigInt> = unsafe {
        ctx.moduli_operators()[0]
            .center_vec_vt(u_rns.coefficients().row(0).as_slice().unwrap())
            .iter()
            .map(|&x| BigInt::from(x))
            .collect()
    };

    let e0: Vec<BigInt> = unsafe {
        ctx.moduli_operators()[0]
            .center_vec_vt(e0_rns.coefficients().row(0).as_slice().unwrap())
            .iter()
            .map(|&x| BigInt::from(x))
            .collect()
    };

    let e1: Vec<BigInt> = unsafe {
        ctx.moduli_operators()[0]
            .center_vec_vt(e1_rns.coefficients().row(0).as_slice().unwrap())
            .iter()
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
    let mut cyclo = vec![BigInt::from(0u64); (N + 1) as usize];
    cyclo[0] = BigInt::from(1u64); // x^0 term
    cyclo[N as usize] = BigInt::from(1u64); // x^(N) term

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
        // Compute the required values for each modulus
        let ct0i: Vec<BigInt> = ct0_coeffs.iter().map(|&x| BigInt::from(x)).collect();
        let ct1i: Vec<BigInt> = ct1_coeffs.iter().map(|&x| BigInt::from(x)).collect();
        let pk0i: Vec<BigInt> = pk0_coeffs.iter().map(|&x| BigInt::from(x)).collect();
        let pk1i: Vec<BigInt> = pk1_coeffs.iter().map(|&x| BigInt::from(x)).collect();

        // k0qi = -t^{-1} mod qi
        let koqi_u64 = qi.inv(qi.neg(t.modulus())).unwrap();
        let k0qi = BigInt::from(koqi_u64);

        // k1 * k0qi as big int
        let ki = scalar_mul(&k1, &k0qi);

        // Calculate ct0i_hat = pk0 * ui + e0i + ki
        let ct0i_hat = {
            let pk0i_times_u = poly_mul(&pk0i, &u);
            assert_eq!((pk0i_times_u.len() as u64) - 1, 2 * (N - 1));

            let e0_plus_ki = poly_add(&e0, &ki);
            assert_eq!((e0_plus_ki.len() as u64) - 1, N - 1);

            poly_add(&pk0i_times_u, &e0_plus_ki)
        };
        assert_eq!((ct0i_hat.len() as u64) - 1, 2 * (N - 1));

        // Compute numerator = ct0i - ct0i_hat
        let ct0i_minus_ct0i_hat = poly_sub(&ct0i, &ct0i_hat);
        assert_eq!((ct0i_minus_ct0i_hat.len() as u64) - 1, 2 * (N - 1));

        // Compute r2i as the quotient of numerator divided by the
        // cyclotomic polynomial, and reduce/center the resulting
        // coefficients to produce: (ct0i - ct0i_hat) / (x^N + 1) mod Z_qi
        let (mut r2i, _) = poly_div(&ct0i_minus_ct0i_hat, &cyclo);
        assert_eq!((r2i.len() as u64) - 1, N - 2); // Order(r2i) = N - 2
        reduce_and_center_coefficients(&mut r2i, &BigInt::from(qi.modulus()));

        // Calculate r1i using r2i and cyclo
        let r1i_num = {
            let r2i_times_cyclo = poly_mul(&r2i, &cyclo);
            assert_eq!((r2i_times_cyclo.len() as u64) - 1, 2 * (N - 1));

            poly_sub(&ct0i_minus_ct0i_hat, &r2i_times_cyclo)
        };
        assert_eq!((r1i_num.len() as u64) - 1, 2 * (N - 1));
        let (r1i, _) = poly_div(&r1i_num, &[BigInt::from(qi.modulus())]);
        assert_eq!((r1i.len() as u64) - 1, 2 * (N - 1)); // Order(r1i) = 2*(n-1)

        /////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Calculate ct1i_hat = pk1 * ui + e1i
        let ct1i_hat = {
            let pk1_ui = poly_mul(&pk1i, &u);
            poly_add(&pk1_ui, &e1)
        };

        // Compute num = ct1i - ct1i_hat
        let mut ct1i_minus_ct1i_hat = poly_sub(&ct1i, &ct1i_hat);

        // Center coefficients of num with respect to the current modulus
        reduce_and_center_coefficients(&mut ct1i_minus_ct1i_hat, &BigInt::from(moduli[i]));

        // Compute p2i as the quotient of num divided by the cyclotomic polynomial
        let (p2i, _) = poly_div(&ct1i_minus_ct1i_hat, &cyclo.clone());

        // Calculate p1i using p2i and cyclo
        let adjusted_num = {
            let p2i_cyclo = poly_mul(&p2i, &cyclo);
            let ct1i_ct1i_hat = poly_sub(&ct1i, &ct1i_hat);
            poly_sub(&ct1i_ct1i_hat, &p2i_cyclo)
        };
        let (p1i, _) = poly_div(&adjusted_num, &[BigInt::from(moduli[i])]);

        p2is.push(p2i);
        p1is.push(p1i);
        r2is.push(r2i);
        r1is.push(r1i);
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

fn center_array_rows(array: ArrayView2<u64>, q: &[Modulus]) -> Array2<i64> {
    let (num_rows, num_cols) = array.dim();

    // Create a new Array2<i64> with zeros
    let mut result = Array2::<i64>::zeros((num_rows, num_cols));

    // Iterate over each row, center it, and assign the result to the new array
    for (i, row) in array.outer_iter().enumerate() {
        let centered_row = unsafe { q[i].center_vec_vt(row.as_slice().unwrap()) };
        result
            .slice_mut(s![i, ..])
            .assign(&Array1::from(centered_row));
    }

    result
}

fn poly_add(poly1: &Vec<BigInt>, poly2: &Vec<BigInt>) -> Vec<BigInt> {
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

fn poly_sub(poly1: &Vec<BigInt>, poly2: &Vec<BigInt>) -> Vec<BigInt> {
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

fn poly_mul(poly1: &Vec<BigInt>, poly2: &Vec<BigInt>) -> Vec<BigInt> {
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
            remainder[i + j] = &remainder[i + j] - &divisor[j] * &coeff;
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
