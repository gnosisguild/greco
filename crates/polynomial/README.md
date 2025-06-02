# Greco Polynomial Library

A Rust library providing polynomial arithmetic modulo operations to support Greco constants generation for correct ciphertext encryption under public key BFV using zero-knowledge proofs.

## Features

- **Basic Polynomial Arithmetic**: Addition, subtraction, multiplication, and division operations
- **Modular Reduction Operations**: Coefficient reduction and centering modulo prime numbers
- **Cyclotomic Polynomial Operations**: Specialized operations for cyclotomic polynomials (x^N + 1)
- **Range Checking**: Validate polynomial coefficients are within specified bounds
- **Arbitrary Precision**: Uses `BigInt` for cryptographic-grade precision
- **Serde Support**: Optional serialization/deserialization support

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
e3-greco-polynomial = "0.1.0"
```

For serialization support:

```toml
[dependencies]
e3-greco-polynomial = { version = "0.1.0", features = ["serde"] }
```

## Quick Start

```rust
use e3_greco_polynomial::{Polynomial, BigInt};

// Create polynomials
let poly1 = Polynomial::new(vec![BigInt::from(2), BigInt::from(3), BigInt::from(1)]); // 2x² + 3x + 1
let poly2 = Polynomial::new(vec![BigInt::from(1), BigInt::from(2)]);                    // x + 2

// Basic arithmetic
let sum = poly1.add(&poly2);
let product = poly1.mul(&poly2);
let difference = poly1.sub(&poly2);

// Division with quotient and remainder
let (quotient, remainder) = poly1.div(&poly2).unwrap();

// Scalar multiplication
let scaled = poly1.scalar_mul(&BigInt::from(5));

// Polynomial evaluation
let result = poly1.evaluate(&BigInt::from(2)); // Evaluate at x = 2
```

## Mathematical Background

This library is designed specifically for cryptographic applications involving polynomial arithmetic in rings of the form `Z_q[X]/(X^N + 1)`, where:

- `Z_q` is the ring of integers modulo a prime `q`
- `X^N + 1` is a cyclotomic polynomial
- Operations preserve the mathematical structure required for BFV homomorphic encryption

### Coefficient Representation

Polynomials are represented with coefficients in **descending order of degree**:

- `[a₂, a₁, a₀]` represents `a₂x² + a₁x + a₀`

### Modular Reduction

The library provides centered reduction where coefficients are mapped to the symmetric range `[-(q-1)/2, (q-1)/2]` rather than the standard range `[0, q-1]`. This is crucial for cryptographic applications where maintaining the sign and magnitude of coefficients is important.

## Error Handling

The library defines a `PolynomialError` enum for handling various error conditions:

```rust
use e3_greco_polynomial::PolynomialError;

match poly1.div(&poly2) {
    Ok((quotient, remainder)) => { /* handle success */ }
    Err(PolynomialError::DivisionByZero) => { /* handle division by zero */ }
    Err(PolynomialError::InvalidPolynomial(msg)) => { /* handle invalid polynomial */ }
    Err(PolynomialError::ModulusError(msg)) => { /* handle modulus error */ }
}
```

## Performance Considerations

- **Arbitrary Precision**: All operations use `BigInt` for unlimited precision, which may be slower than native integer types but ensures correctness for cryptographic applications
- **Memory Usage**: Large degree polynomials with many coefficients will consume significant memory
- **Algorithm Complexity**:
  - Addition/Subtraction: O(n)
  - Multiplication: O(n²) using naive algorithm
  - Division: O(n²) using polynomial long division

## Use Cases

This library is specifically designed for:

- **Zero-Knowledge Proofs**: Generating polynomial commitments and proofs
- **Homomorphic Encryption**: BFV scheme parameter generation
- **Lattice Cryptography**: Polynomial operations in cyclotomic rings
- **Cryptographic Research**: Experimenting with polynomial-based cryptographic primitives

## Contributing

Contributions are welcome! Please ensure that:

1. All tests pass: `cargo test`
2. Code is properly formatted: `cargo fmt`
3. No clippy warnings: `cargo clippy`
4. New features include comprehensive tests
5. Documentation is updated for API changes

## License

This project is licensed under the MIT License - see the [LICENSE](../../LICENSE) file for details.
