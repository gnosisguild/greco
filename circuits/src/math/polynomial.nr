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
        let mut acc = self.coefficients[0];
        for i in 1..self.coefficients.len() {
            acc = acc * x + self.coefficients[i];
        }
        acc
    }

    ///  Adds `upper_bound` to the coefficients of the polynomial and constrains them to be in the range `[0, 2 * upper_bound]`.
    pub fn range_check_1bound(self, upper_bound: u64) {
        let comp_value = (2 * upper_bound) as Field;
        for i in 0..self.coefficients.len() {
            let shifted_coeff = self.coefficients[i] + upper_bound as Field;
            shifted_coeff.assert_max_bit_size::<66>();
            (comp_value - shifted_coeff).assert_max_bit_size::<66>();
        }
    }

    /// Adds `-lower_bound` to the coefficients of the polynomial and constrains them to be in the range `[0, upper_bound - lower_bound]`.
    pub fn range_check_2bounds(self, upper_bound: u64, lower_bound: i64) {
        let mut value_shift = (lower_bound * -1) as u64 as Field;
        let comp_value = value_shift + upper_bound as Field;
        for i in 0..self.coefficients.len() {
            let shifted_coeff = self.coefficients[i] + value_shift;
            shifted_coeff.assert_max_bit_size::<66>();
            (comp_value - shifted_coeff).assert_max_bit_size::<66>();
        }
    }
}

#[test]
fn test_polynomial_eval() {
    let coeffs = [1, 2, 3]; // represents 1x^2 + 2x + 3
    let poly = Polynomial::new(coeffs);

    let x = 2; // evaluate at x = 2
    let result = poly.eval(x);
    // (1 * 2^2) + (2 * 2) + 3 = 4 + 4 + 3 = 11
    assert_eq(result, 11);
}

#[test]
fn test_polynomial_bounds() {
    let coeffs = [-16, 240, 242];
    let poly = Polynomial::new(coeffs);

    // Test single bound check
    poly.range_check_1bound(242);

    // Test double bounds check
    poly.range_check_2bounds(242, -16);
}

#[test]
fn test_polynomial_eval_zero() {
    let coeffs = [1, -2, 1]; // x^2 - 2x + 1 = (x-1)^2
    let poly = Polynomial::new(coeffs);

    let x = 1; // evaluate at x = 1, should be 0
    let result = poly.eval(x);
    assert_eq(result, 0);
}

#[test]
fn test_polynomial_eval_large() {
    let coeffs = [1000000, -500000, 250000]; // large coefficients
    let poly = Polynomial::new(coeffs);

    let x = 2;
    let result = poly.eval(x);
    // 1000000*2^2 - 500000*2 + 250000 = 4000000 - 1000000 + 250000 = 3250000
    assert_eq(result, 3250000);
}
