use std::field::bn254::assert_lt;

/// Polynomial structure
pub struct Polynomial<let N: u32> {
    // Coefficients of the polynomial (Biggest degree at the first index)
    pub coefficients: [Field; N],
}

impl<let N: u32> Polynomial<N> {
    /// Creates a new polynomial
    pub fn new(coefficients: [Field; N]) -> Self {
        Polynomial { coefficients }
    }

    /// Evaluate a polynomial at a point in the field.
    pub fn eval(self, x: Field) -> Field {
        {
            let mut acc = self.coefficients[0];
            for i in 1..self.coefficients.len() {
                acc = acc * x + self.coefficients[i];
            }
            acc
        }
    }

    ///  Adds `upper_bound` to the coefficients of the polynomial and constrains them to be in the range `[0, 2 * upper_bound]`.
    pub fn range_check_1bound(self, upper_bound: u64) {
        for i in 0..self.coefficients.len() {
            let shifted_coeff = self.coefficients[i] + upper_bound as Field;
            assert_lt(shifted_coeff, (2 * upper_bound as Field) + 1);
        }
    }

    /// Adds `-lower_bound` to the coefficients of the polynomial and constrains them to be in the range `[0, upper_bound - lower_bound]`.
    pub fn range_check_2bounds(self, upper_bound: u64, lower_bound: i64) {
        let mut value_shift = (lower_bound * -1) as u64;
        for i in 0..self.coefficients.len() {
            let shifted_coeff = self.coefficients[i] + value_shift as Field;
            assert_lt(
                shifted_coeff,
                upper_bound as Field - lower_bound as Field + 1,
            );
        }
    }
}

#[test]
// Tests single bound check and double bounds check with a small Polynomial
fn test_polynomial_bounds() {
    let mut coeffs = [-16, 240, 242];

    let mut poly = Polynomial::new(coeffs);
    poly.range_check_1bound(242);
    poly.range_check_2bounds(242, -16);
}

#[test]
// Tests polynomial evaluation on a random small number
fn test_polynomial_eval() {
    let mut coeffs = [1, 2, 3];
    let mut poly = Polynomial::new(coeffs);

    let eval = poly.eval(3);
    assert_eq(eval as u64, 18);
}
