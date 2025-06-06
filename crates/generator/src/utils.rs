//! Utility functions for the Greco generator.
//!
//! This module contains helper functions for string conversion,
//! serialization, and other common operations.

use num_bigint::BigInt;

/// Convert a 1D vector of BigInt to a vector of strings.
pub fn to_string_1d_vec(poly: &Vec<BigInt>) -> Vec<String> {
    poly.iter().map(|coef| coef.to_string()).collect()
}

/// Convert a 2D vector of BigInt to a vector of vectors of strings.
pub fn to_string_2d_vec(poly: &Vec<Vec<BigInt>>) -> Vec<Vec<String>> {
    poly.iter().map(|row| to_string_1d_vec(row)).collect()
}
