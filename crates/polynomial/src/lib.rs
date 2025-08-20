//! # E3 Greco Polynomial Library
//!
//! This library provides polynomial arithmetic modulo operations
//! to support Greco constants generation for correct ciphertext encryption
//! under public key BFV using zero-knowledge proofs.
//!
//! ## Features
//!
//! - Basic polynomial arithmetic (add, subtract, multiply, divide)
//! - Modular reduction operations
//! - Cyclotomic polynomial operations
//! - Coefficient centering and range checking
//!
//! ## Example
//!
//! ```rust
//! use e3_greco_polynomial::{Polynomial, BigInt};
//!
//! let poly1 = Polynomial::new(vec![BigInt::from(1), BigInt::from(2), BigInt::from(3)]);
//! let poly2 = Polynomial::new(vec![BigInt::from(1), BigInt::from(1)]);
//! let result = poly1.add(&poly2);
//! ```

/// Re-export for convenience
pub use num_bigint::BigInt;
pub use num_traits::{One, Zero};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A polynomial represented by its coefficients in descending order of degree.
///
/// The coefficients are stored as `BigInt` to support arbitrary precision arithmetic
/// required for cryptographic operations.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Polynomial {
    /// Coefficients in descending order (highest degree first)
    coefficients: Vec<BigInt>,
}

/// Errors that can occur during polynomial operations
#[derive(Debug, Clone, PartialEq)]
pub enum PolynomialError {
    /// Division by zero polynomial
    DivisionByZero,
    /// Invalid polynomial (e.g., empty coefficients)
    InvalidPolynomial(String),
    /// Modulus operation error
    ModulusError(String),
}

impl std::fmt::Display for PolynomialError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PolynomialError::DivisionByZero => write!(f, "Division by zero polynomial"),
            PolynomialError::InvalidPolynomial(msg) => write!(f, "Invalid polynomial: {}", msg),
            PolynomialError::ModulusError(msg) => write!(f, "Modulus error: {}", msg),
        }
    }
}

impl std::error::Error for PolynomialError {}

impl Polynomial {
    /// Creates a new polynomial from a vector of coefficients.
    ///
    /// # Arguments
    ///
    /// * `coefficients` - Vector of coefficients in descending order of degree
    ///
    /// # Examples
    ///
    /// ```rust
    /// use e3_greco_polynomial::{Polynomial, BigInt};
    ///
    /// // Creates polynomial: 2x^2 + 3x + 1
    /// let poly = Polynomial::new(vec![BigInt::from(2), BigInt::from(3), BigInt::from(1)]);
    /// ```
    pub fn new(coefficients: Vec<BigInt>) -> Self {
        Self { coefficients }
    }

    /// Creates a zero polynomial of specified degree.
    ///
    /// # Arguments
    ///
    /// * `degree` - The degree of the zero polynomial
    ///
    /// # Examples
    ///
    /// ```rust
    /// use e3_greco_polynomial::Polynomial;
    ///
    /// let zero_poly = Polynomial::zero(3);
    /// ```
    pub fn zero(degree: usize) -> Self {
        Self {
            coefficients: vec![BigInt::zero(); degree + 1],
        }
    }

    /// Creates a constant polynomial.
    ///
    /// # Arguments
    ///
    /// * `constant` - The constant value
    ///
    /// # Examples
    ///
    /// ```rust
    /// use e3_greco_polynomial::{Polynomial, BigInt};
    ///
    /// let const_poly = Polynomial::constant(BigInt::from(42));
    /// ```
    pub fn constant(constant: BigInt) -> Self {
        Self {
            coefficients: vec![constant],
        }
    }

    /// Returns the coefficients of the polynomial.
    pub fn coefficients(&self) -> &[BigInt] {
        &self.coefficients
    }

    /// Returns the degree of the polynomial.
    pub fn degree(&self) -> usize {
        if self.coefficients.is_empty() {
            0
        } else {
            self.coefficients.len() - 1
        }
    }

    /// Checks if the polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.coefficients.iter().all(|c| c.is_zero())
    }

    /// Removes leading zero coefficients from the polynomial.
    pub fn trim_leading_zeros(mut self) -> Self {
        while self.coefficients.len() > 1 && self.coefficients[0].is_zero() {
            self.coefficients.remove(0);
        }
        self
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
    /// * `other` - A reference to the polynomial to add to `self`
    ///
    /// # Returns
    ///
    /// A new polynomial containing the sum of the two polynomials
    ///
    /// # Examples
    ///
    /// ```rust
    /// use e3_greco_polynomial::{Polynomial, BigInt};
    ///
    /// let poly1 = Polynomial::new(vec![BigInt::from(1), BigInt::from(2)]);
    /// let poly2 = Polynomial::new(vec![BigInt::from(3), BigInt::from(4)]);
    /// let result = poly1.add(&poly2);
    /// ```
    pub fn add(&self, other: &Self) -> Self {
        let max_length = std::cmp::max(self.coefficients.len(), other.coefficients.len());
        let mut result = vec![BigInt::zero(); max_length];

        // Copy coefficients from the first polynomial
        for (i, coeff) in self.coefficients.iter().enumerate() {
            result[max_length - self.coefficients.len() + i] = coeff.clone();
        }

        // Add coefficients from the second polynomial
        for (i, coeff) in other.coefficients.iter().enumerate() {
            result[max_length - other.coefficients.len() + i] += coeff;
        }

        Polynomial::new(result)
    }

    /// Subtracts one polynomial from another.
    ///
    /// This function subtracts the second polynomial from the first polynomial by
    /// negating the coefficients of the second polynomial and then adding.
    ///
    /// # Arguments
    ///
    /// * `other` - A reference to the polynomial to subtract from `self`
    ///
    /// # Returns
    ///
    /// A new polynomial containing the difference
    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    /// Negates all coefficients of the polynomial.
    ///
    /// # Returns
    ///
    /// A new polynomial with all coefficients negated
    pub fn neg(&self) -> Self {
        Polynomial::new(self.coefficients.iter().map(|x| -x).collect())
    }

    /// Multiplies two polynomials using naive algorithm.
    ///
    /// Given two polynomials, this function computes their product using the
    /// standard convolution algorithm.
    ///
    /// # Arguments
    ///
    /// * `other` - A reference to the polynomial to multiply with `self`
    ///
    /// # Returns
    ///
    /// A new polynomial containing the product
    pub fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Polynomial::zero(0);
        }

        let product_len = self.coefficients.len() + other.coefficients.len() - 1;
        let mut product = vec![BigInt::zero(); product_len];

        for i in 0..self.coefficients.len() {
            for j in 0..other.coefficients.len() {
                product[i + j] += &self.coefficients[i] * &other.coefficients[j];
            }
        }

        Polynomial::new(product)
    }

    /// Divides one polynomial by another, returning the quotient and remainder.
    ///
    /// This function performs polynomial long division.
    ///
    /// # Arguments
    ///
    /// * `divisor` - A reference to the divisor polynomial
    ///
    /// # Returns
    ///
    /// A result containing a tuple of (quotient, remainder) or an error
    ///
    /// # Errors
    ///
    /// Returns `PolynomialError::DivisionByZero` if the divisor is zero
    /// Returns `PolynomialError::InvalidPolynomial` if the divisor has zero leading coefficient
    pub fn div(&self, divisor: &Self) -> Result<(Self, Self), PolynomialError> {
        if divisor.is_zero() {
            return Err(PolynomialError::DivisionByZero);
        }

        if divisor.coefficients.is_empty() || divisor.coefficients[0].is_zero() {
            return Err(PolynomialError::InvalidPolynomial(
                "Leading coefficient of divisor cannot be zero".to_string(),
            ));
        }

        if self.degree() < divisor.degree() {
            return Ok((Polynomial::zero(0), self.clone()));
        }

        let mut quotient =
            vec![BigInt::zero(); self.coefficients.len() - divisor.coefficients.len() + 1];
        let mut remainder = self.coefficients.clone();

        for i in 0..quotient.len() {
            if i >= remainder.len() {
                break;
            }
            let coeff = &remainder[i] / &divisor.coefficients[0];
            quotient[i] = coeff.clone();

            for j in 0..divisor.coefficients.len() {
                if i + j < remainder.len() {
                    remainder[i + j] = &remainder[i + j] - &divisor.coefficients[j] * &coeff;
                }
            }
        }

        // Remove leading zero coefficients from remainder
        while !remainder.is_empty() && remainder[0].is_zero() {
            remainder.remove(0);
        }

        Ok((Polynomial::new(quotient), Polynomial::new(remainder)))
    }

    /// Multiplies each coefficient of the polynomial by a scalar.
    ///
    /// # Arguments
    ///
    /// * `scalar` - A `BigInt` scalar to multiply with each coefficient
    ///
    /// # Returns
    ///
    /// A new polynomial with each coefficient multiplied by the scalar
    pub fn scalar_mul(&self, scalar: &BigInt) -> Self {
        Polynomial::new(self.coefficients.iter().map(|x| x * scalar).collect())
    }

    /// Reduces the polynomial modulo a cyclotomic polynomial.
    ///
    /// This function performs polynomial division by the cyclotomic polynomial
    /// and returns the remainder.
    ///
    /// # Arguments
    ///
    /// * `cyclo` - Coefficients of the cyclotomic polynomial
    ///
    /// # Returns
    ///
    /// A new polynomial representing the remainder after reduction
    pub fn reduce_by_cyclotomic(&self, cyclo: &[BigInt]) -> Result<Self, PolynomialError> {
        let cyclo_poly = Polynomial::new(cyclo.to_vec());
        let (_, remainder) = self.div(&cyclo_poly)?;

        let n = cyclo.len() - 1;
        let mut out = vec![BigInt::zero(); n];

        if !remainder.coefficients.is_empty() {
            let start_idx = n.saturating_sub(remainder.coefficients.len());
            let end_idx = std::cmp::min(start_idx + remainder.coefficients.len(), n);
            out[start_idx..end_idx]
                .clone_from_slice(&remainder.coefficients[..end_idx - start_idx]);
        }

        Ok(Polynomial::new(out))
    }

    /// Reduces coefficients modulo a prime and centers them.
    ///
    /// # Arguments
    ///
    /// * `modulus` - The prime modulus
    ///
    /// # Returns
    ///
    /// A new polynomial with coefficients reduced and centered
    pub fn reduce_and_center(&self, modulus: &BigInt) -> Self {
        let half_modulus = modulus / 2;
        let reduced_coeffs = self
            .coefficients
            .iter()
            .map(|x| reduce_and_center(x, modulus, &half_modulus))
            .collect();
        Polynomial::new(reduced_coeffs)
    }

    /// Evaluates the polynomial at a given point.
    ///
    /// # Arguments
    ///
    /// * `x` - The point at which to evaluate the polynomial
    ///
    /// # Returns
    ///
    /// The value of the polynomial at the given point
    pub fn evaluate(&self, x: &BigInt) -> BigInt {
        if self.coefficients.is_empty() {
            return BigInt::zero();
        }

        // Use Horner's method for efficient evaluation
        let mut result = self.coefficients[0].clone();
        for coeff in &self.coefficients[1..] {
            result = result * x + coeff;
        }
        result
    }
}

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

/// Reduces and centers a scalar 
pub fn reduce_and_center_scalar(x: &BigInt, modulus: &BigInt) -> BigInt {
    let half_modulus = modulus / 2;
    reduce_and_center(x, modulus, &half_modulus)
}
///
/// Reduces a scalar 
pub fn reduce_scalar(x: &BigInt, modulus: &BigInt) -> BigInt {
    let r = (x + modulus) % modulus;
    r
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
    let reduced = poly
        .reduce_by_cyclotomic(cyclo)
        .expect("Failed to reduce by cyclotomic");
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

/// Checks if all coefficients in a vector are within a centered range.
///
/// This function verifies that every coefficient in the input vector falls within
/// the inclusive range [lower_bound, upper_bound]. This is typically used for
/// coefficients that have been centered around zero.
///
/// # Arguments
///
/// * `vec` - A slice of `BigInt` coefficients to check
/// * `lower_bound` - The minimum allowed value (inclusive)
/// * `upper_bound` - The maximum allowed value (inclusive)
///
/// # Returns
///
/// * `true` if all coefficients are within bounds, `false` otherwise
///
/// # Examples
///
/// ```
/// use e3_greco_polynomial::{range_check_centered, BigInt};
///
/// let coeffs = vec![BigInt::from(-2), BigInt::from(0), BigInt::from(2)];
/// let result = range_check_centered(&coeffs, &BigInt::from(-3), &BigInt::from(3));
/// assert!(result);
/// ```
pub fn range_check_centered(vec: &[BigInt], lower_bound: &BigInt, upper_bound: &BigInt) -> bool {
    vec.iter()
        .all(|coeff| coeff >= lower_bound && coeff <= upper_bound)
}

/// Checks if all coefficients satisfy standard range constraints with separate upper and lower bounds.
///
/// This function verifies that each coefficient falls within one of two ranges:
/// 1. [0, up_bound] (positive range)
/// 2. [modulus + low_bound, modulus) (negative range wrapped around modulus)
///
/// This is commonly used in cryptographic applications where coefficients can be
/// represented in both positive and negative forms modulo a prime.
///
/// # Arguments
///
/// * `vec` - A slice of `BigInt` coefficients to check
/// * `low_bound` - The lower bound for the negative range (typically negative)
/// * `up_bound` - The upper bound for the positive range
/// * `modulus` - The modulus used for wraparound calculations
///
/// # Returns
///
/// * `true` if all coefficients satisfy the range constraints, `false` otherwise
///
/// # Mathematical Background
///
/// In modular arithmetic, negative values are often represented as their positive
/// equivalents: `-x ≡ modulus - x (mod modulus)`. This function checks both
/// the direct positive representation and the wrapped negative representation.
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

/// Checks if all coefficients satisfy symmetric standard range constraints.
///
/// This function verifies that each coefficient falls within one of two symmetric ranges:
/// 1. [0, bound] (positive range)
/// 2. [modulus - bound, modulus) (negative range wrapped around modulus)
///
/// This is a special case of `range_check_standard_2bounds` where the bounds are
/// symmetric around zero. Commonly used for error distributions in cryptography.
///
/// # Arguments
///
/// * `vec` - A slice of `BigInt` coefficients to check
/// * `bound` - The symmetric bound (both positive and negative)
/// * `modulus` - The modulus used for wraparound calculations
///
/// # Returns
///
/// * `true` if all coefficients satisfy the symmetric range constraints, `false` otherwise
///
/// # Examples
///
/// ```
/// use e3_greco_polynomial::{range_check_standard, BigInt};
///
/// let coeffs = vec![BigInt::from(3), BigInt::from(0)];
/// let result = range_check_standard(&coeffs, &BigInt::from(5), &BigInt::from(7));
/// assert!(result);
/// ```
///
/// # Mathematical Background
///
/// For a coefficient `c` and bound `b`, this function accepts:
/// - `c ∈ [0, b]` (small positive values)
/// - `c ∈ [modulus - b, modulus)` (small negative values as positive representatives)
pub fn range_check_standard(vec: &[BigInt], bound: &BigInt, modulus: &BigInt) -> bool {
    vec.iter().all(|coeff| {
        (coeff >= &BigInt::from(0) && coeff <= bound)
            || (coeff >= &(modulus - bound) && coeff < modulus)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polynomial_creation() {
        let poly = Polynomial::new(vec![BigInt::from(1), BigInt::from(2), BigInt::from(3)]);
        assert_eq!(
            poly.coefficients(),
            &[BigInt::from(1), BigInt::from(2), BigInt::from(3)]
        );
        assert_eq!(poly.degree(), 2);
    }

    #[test]
    fn test_zero_polynomial() {
        let zero = Polynomial::zero(3);
        assert_eq!(zero.degree(), 3);
        assert!(zero.is_zero());
    }

    #[test]
    fn test_constant_polynomial() {
        let const_poly = Polynomial::constant(BigInt::from(42));
        assert_eq!(const_poly.degree(), 0);
        assert_eq!(const_poly.coefficients(), &[BigInt::from(42)]);
    }

    #[test]
    fn test_polynomial_addition() {
        let poly1 = Polynomial::new(vec![BigInt::from(1), BigInt::from(2)]);
        let poly2 = Polynomial::new(vec![BigInt::from(3), BigInt::from(4)]);
        let result = poly1.add(&poly2);
        assert_eq!(result.coefficients(), &[BigInt::from(4), BigInt::from(6)]);
    }

    #[test]
    fn test_polynomial_subtraction() {
        let poly1 = Polynomial::new(vec![BigInt::from(5), BigInt::from(3)]);
        let poly2 = Polynomial::new(vec![BigInt::from(2), BigInt::from(1)]);
        let result = poly1.sub(&poly2);
        assert_eq!(result.coefficients(), &[BigInt::from(3), BigInt::from(2)]);
    }

    #[test]
    fn test_polynomial_negation() {
        let poly = Polynomial::new(vec![BigInt::from(1), BigInt::from(-2), BigInt::from(3)]);
        let neg_poly = poly.neg();
        assert_eq!(
            neg_poly.coefficients(),
            &[BigInt::from(-1), BigInt::from(2), BigInt::from(-3)]
        );
    }

    #[test]
    fn test_polynomial_multiplication() {
        let poly1 = Polynomial::new(vec![BigInt::from(1), BigInt::from(2)]); // x + 2
        let poly2 = Polynomial::new(vec![BigInt::from(1), BigInt::from(3)]); // x + 3
        let result = poly1.mul(&poly2); // Should be x^2 + 5x + 6
        assert_eq!(
            result.coefficients(),
            &[BigInt::from(1), BigInt::from(5), BigInt::from(6)]
        );
    }

    #[test]
    fn test_polynomial_division() {
        let dividend = Polynomial::new(vec![BigInt::from(1), BigInt::from(5), BigInt::from(6)]); // x^2 + 5x + 6
        let divisor = Polynomial::new(vec![BigInt::from(1), BigInt::from(2)]); // x + 2
        let (quotient, remainder) = dividend.div(&divisor).unwrap();
        assert_eq!(quotient.coefficients(), &[BigInt::from(1), BigInt::from(3)]); // x + 3
        assert!(remainder.is_zero());
    }

    #[test]
    fn test_division_by_zero() {
        let poly = Polynomial::new(vec![BigInt::from(1), BigInt::from(2)]);
        let zero = Polynomial::zero(0);
        assert!(matches!(
            poly.div(&zero),
            Err(PolynomialError::DivisionByZero)
        ));
    }

    #[test]
    fn test_scalar_multiplication() {
        let poly = Polynomial::new(vec![BigInt::from(1), BigInt::from(2), BigInt::from(3)]);
        let scalar = BigInt::from(5);
        let result = poly.scalar_mul(&scalar);
        assert_eq!(
            result.coefficients(),
            &[BigInt::from(5), BigInt::from(10), BigInt::from(15)]
        );
    }

    #[test]
    fn test_polynomial_evaluation() {
        let poly = Polynomial::new(vec![BigInt::from(1), BigInt::from(2), BigInt::from(3)]); // x^2 + 2x + 3
        let result = poly.evaluate(&BigInt::from(2)); // 1*4 + 2*2 + 3 = 11
        assert_eq!(result, BigInt::from(11));
    }

    #[test]
    fn test_trim_leading_zeros() {
        let poly = Polynomial::new(vec![
            BigInt::from(0),
            BigInt::from(0),
            BigInt::from(1),
            BigInt::from(2),
        ]);
        let trimmed = poly.trim_leading_zeros();
        assert_eq!(trimmed.coefficients(), &[BigInt::from(1), BigInt::from(2)]);
    }

    #[test]
    fn test_reduce_and_center() {
        let poly = Polynomial::new(vec![BigInt::from(10), BigInt::from(15), BigInt::from(20)]);
        let modulus = BigInt::from(7);
        let result = poly.reduce_and_center(&modulus);
        // 10 % 7 = 3, 15 % 7 = 1, 20 % 7 = 6 -> -1 (centered)
        assert_eq!(
            result.coefficients(),
            &[BigInt::from(3), BigInt::from(1), BigInt::from(-1)]
        );
    }

    #[test]
    fn test_reduce_and_center_function() {
        let modulus = BigInt::from(7);
        let half_modulus = &modulus / 2;

        // Test positive number
        assert_eq!(
            reduce_and_center(&BigInt::from(10), &modulus, &half_modulus),
            BigInt::from(3)
        );

        // Test negative number
        assert_eq!(
            reduce_and_center(&BigInt::from(-3), &modulus, &half_modulus),
            BigInt::from(-3)
        );

        // Test number greater than half modulus
        assert_eq!(
            reduce_and_center(&BigInt::from(6), &modulus, &half_modulus),
            BigInt::from(-1)
        );
    }

    #[test]
    fn test_reduce_coefficients() {
        let coeffs = vec![BigInt::from(10), BigInt::from(-3), BigInt::from(15)];
        let modulus = BigInt::from(7);
        let result = reduce_coefficients(&coeffs, &modulus);
        assert_eq!(
            result,
            vec![BigInt::from(3), BigInt::from(4), BigInt::from(1)]
        );
    }

    #[test]
    fn test_range_check_centered() {
        let vec = vec![BigInt::from(-2), BigInt::from(0), BigInt::from(2)];
        let lower = BigInt::from(-3);
        let upper = BigInt::from(3);
        assert!(range_check_centered(&vec, &lower, &upper));

        let vec_out_of_range = vec![BigInt::from(-5), BigInt::from(0), BigInt::from(2)];
        assert!(!range_check_centered(&vec_out_of_range, &lower, &upper));
    }

    #[test]
    fn test_range_check_standard() {
        let vec = vec![BigInt::from(1), BigInt::from(2), BigInt::from(3)];
        let bound = BigInt::from(5);
        let modulus = BigInt::from(7);
        assert!(range_check_standard(&vec, &bound, &modulus));
    }
}
