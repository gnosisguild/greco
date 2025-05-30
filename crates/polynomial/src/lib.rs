/// @todo needs cleanup and refactoring and tests.

/// Provides helper methods that perform modular poynomial arithmetic over polynomials encoded in vectors
/// of coefficients from largest degree to lowest.
use num_bigint::BigInt;
use num_traits::*;

struct Polynomial {
    coefficients: Vec<BigInt>,
}

/// @todo cont. still tbd
impl Polynomial {
    pub fn new(coefficients: Vec<BigInt>) -> Self {
        Self { coefficients }
    }

    /// Adds two polynomials together.
    ///
    /// This function performs polynomial addition by:
    /// 1. Finding the maximum length between the two polynomials
    /// 2. Creating a new polynomial with the maximum length
    /// 3. Adding the coefficients of both polynomials term by term
    ///
    /// # Arguments
    ///
    /// * `p` - A reference to the polynomial to add to `self`
    ///
    /// # Returns
    ///
    /// A new polynomial containing the sum of the two polynomials
    pub fn add(&self, p: &Self) -> Self {
        let max_length = std::cmp::max(self.coefficients.len(), p.coefficients.len());
        let mut result = vec![BigInt::zero(); max_length];

        // Copy coefficients from the first polynomial
        for (i, coeff) in self.coefficients.iter().enumerate() {
            result[max_length - self.coefficients.len() + i] = coeff.clone();
        }

        // Add coefficients from the second polynomial
        for (i, coeff) in p.coefficients.iter().enumerate() {
            result[max_length - p.coefficients.len() + i] += coeff;
        }

        Polynomial::new(result)
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
    pub fn sub(&self, p: &Self) -> Self {
        self.add(&p.neg())
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
    pub fn neg(&self) -> Self {
        Polynomial::new(self.coefficients.iter().map(|x| -x).collect())
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
    pub fn mul(&self, p: &Self) -> Self {
        let product_len = self.coefficients.len() + p.coefficients.len() - 1;
        let mut product = vec![BigInt::zero(); product_len];

        for i in 0..self.coefficients.len() {
            for j in 0..p.coefficients.len() {
                product[i + j] += &self.coefficients[i] * &p.coefficients[j];
            }
        }

        Polynomial::new(product)
    }

    /// Divides one polynomial by another, returning the quotient and remainder.
    ///
    /// This function performs polynomial long division by:
    /// 1. Checking that the divisor is valid (non-empty and leading coefficient non-zero)
    /// 2. Computing the quotient and remainder using the standard long division algorithm
    ///
    /// # Arguments
    ///
    /// * `p` - A reference to the divisor polynomial
    ///
    /// # Returns
    ///
    /// A tuple containing:
    /// * The quotient polynomial
    /// * The remainder polynomial
    ///
    /// # Panics
    ///
    /// This function will panic if:
    /// * The divisor is empty
    /// * The leading coefficient of the divisor is zero
    pub fn div(&self, p: &Self) -> (Self, Self) {
        assert!(
            !p.coefficients.is_empty() && !p.coefficients[0].is_zero(),
            "Leading coefficient of divisor cannot be zero"
        );

        let mut quotient = vec![BigInt::zero(); self.coefficients.len() - p.coefficients.len() + 1];
        let mut remainder = self.coefficients.clone();

        for i in 0..quotient.len() {
            let coeff = &remainder[i] / &p.coefficients[0];
            quotient[i] = coeff.clone();

            for j in 0..p.coefficients.len() {
                remainder[i + j] = &remainder[i + j] - &p.coefficients[j] * &coeff;
            }
        }

        // Remove leading zero coefficients from remainder
        while !remainder.is_empty() && remainder[0].is_zero() {
            remainder.remove(0);
        }

        (Polynomial::new(quotient), Polynomial::new(remainder))
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
    pub fn scalar_mul(&self, scalar: &BigInt) -> Self {
        Polynomial::new(self.coefficients.iter().map(|x| x * scalar).collect())
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
    pub fn reduce_coefficients_by_cyclo(&self, cyclo: &[BigInt]) -> Self {
        let (_, remainder) = self.div(&Polynomial::new(cyclo.to_vec()));

        let N = cyclo.len() - 1;
        let mut out: Vec<BigInt> = vec![BigInt::zero(); N];
        let start_idx = N - remainder.coefficients.len();
        out[start_idx..].clone_from_slice(&remainder.coefficients);

        Polynomial::new(out)
    }
}

/// @todo cont.

/// Reduces a number modulo a prime modulus and centers it.
///
/// This function takes an arbitrary number and reduces it modulo the specified prime modulus.
/// After reduction, the number is adjusted to be within the symmetric range
/// [(−(modulus−1))/2, (modulus−1)/2]. If the number is already within this range, it remains unchanged.
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
pub fn reduce_and_center(x: &BigInt, modulus: &BigInt, half_modulus: &BigInt) -> BigInt {
    // Calculate the remainder ensuring it's non-negative
    let mut r: BigInt = x % modulus;
    if r < BigInt::zero() {
        r += modulus;
    }

    // Adjust the remainder if it is greater than half_modulus
    if (modulus % BigInt::from(2)) == BigInt::from(1) {
        if r > *half_modulus {
            r -= modulus;
        }
    } else if r >= *half_modulus {
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
pub fn reduce_and_center_coefficients_mut(coefficients: &mut [BigInt], modulus: &BigInt) {
    let half_modulus = modulus / BigInt::from(2);
    coefficients
        .iter_mut()
        .for_each(|x| *x = reduce_and_center(x, modulus, &half_modulus));
}

pub fn reduce_and_center_coefficients(
    coefficients: &mut [BigInt],
    modulus: &BigInt,
) -> Vec<BigInt> {
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
pub fn reduce_in_ring(coefficients: &mut Vec<BigInt>, cyclo: &[BigInt], modulus: &BigInt) {
    let poly = Polynomial::new(coefficients.clone());
    let reduced = poly.reduce_coefficients_by_cyclo(cyclo);
    *coefficients = reduced.coefficients;
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
pub fn reduce_coefficients(coefficients: &[BigInt], p: &BigInt) -> Vec<BigInt> {
    coefficients.iter().map(|coeff| (coeff + p) % p).collect()
}

pub fn reduce_coefficients_2d(coefficient_matrix: &[Vec<BigInt>], p: &BigInt) -> Vec<Vec<BigInt>> {
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
pub fn reduce_coefficients_mut(coefficients: &mut [BigInt], p: &BigInt) {
    for coeff in coefficients.iter_mut() {
        *coeff += p;
        *coeff %= p;
    }
}

pub fn range_check_centered(vec: &[BigInt], lower_bound: &BigInt, upper_bound: &BigInt) -> bool {
    vec.iter()
        .all(|coeff| coeff >= lower_bound && coeff <= upper_bound)
}

pub fn range_check_standard_2bounds(
    vec: &[BigInt],
    low_bound: &BigInt,
    up_bound: &BigInt,
    modulus: &BigInt,
) -> bool {
    vec.iter().all(|coeff| {
        (coeff >= &BigInt::from(0) && coeff <= up_bound)
            || (coeff >= &(modulus + low_bound) && coeff < modulus)
    })
}

pub fn range_check_standard(vec: &[BigInt], bound: &BigInt, modulus: &BigInt) -> bool {
    vec.iter().all(|coeff| {
        (coeff >= &BigInt::from(0) && coeff <= bound)
            || (coeff >= &(modulus - bound) && coeff < modulus)
    })
}
