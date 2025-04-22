use std::field::bn254::assert_lt;

/// Polynomial structure
pub struct Polynomial {
    // Coefficients of the polynomial (Biggest degree at the first index)
    pub coefficients: Vec<Field>,
}

impl Polynomial {
    /// Creates a new polynomial
    pub fn new(coefficients: Vec<Field>) -> Self {
        Polynomial { coefficients }
    }

    /// Evaluate a polynomial at a point in the field.
    pub fn eval(self, x: Field) -> Field {
        {
            let mut acc = self.coefficients.get(0);
            for i in 1..self.coefficients.len() {
                acc = acc * x + self.coefficients.get(i);
            }
            acc
        }
    }

    /// Adds `upper_bound` to the coefficients of the polynomial and constrains them to be in the range `[0, 2*upper_bound]`.
    pub fn range_check_1bound(self, upper_bound: u64) {
        for coeff in 0..self.coefficients.len() {
            assert_lt(self.coefficients.get(coeff), 2 * upper_bound as Field + 1);
        }
    }

    /// Adds `-lower_bound` to the coefficients of the polynomial and constrains them to be in the range `[0, upper_bound - lower_bound]`.
    pub fn range_check_2bounds(self, upper_bound: u64, lower_bound: i64) {
        let mut value_shift = lower_bound as u64;
        if lower_bound < 0 {
            value_shift = (lower_bound * -1) as u64;
        }
        for coeff in 0..self.coefficients.len() {
            assert_lt(
                self.coefficients.get(coeff),
                upper_bound as Field - lower_bound as Field + 1,
            );
        }
    }
}

#[test]
fn test_polynomial_bounds() {
    let mut coeffs = Vec::new();
    coeffs.push(23);
    coeffs.push(484);
    coeffs.push(10);

    let mut poly = Polynomial::new(coeffs);

    poly.range_check_1bound(242);
    poly.range_check_2bounds(500, -16);
}

#[test]
fn test_polynomial_eval() {
    let mut coeffs = Vec::new();
    coeffs.push(1);
    coeffs.push(2);
    coeffs.push(3);
    let mut poly = Polynomial::new(coeffs);

    let eval = poly.eval(3);
    assert_eq(eval as u64, 18);
}
