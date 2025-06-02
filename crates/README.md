# Greco Rust Crates

This directory contains the refactored Rust crates for the Greco project, providing modular libraries for polynomial arithmetic and cryptographic parameter generation for BFV encryption proofs.

## Overview

The crates are organized into two main libraries:

- **`e3-greco-polynomial`**: Core polynomial arithmetic library
- **`e3-greco-generator`**: Parameter generation and validation library with CLI

## Crates

### 📐 e3-greco-polynomial

A comprehensive polynomial arithmetic library supporting:

- **Basic Operations**: Addition, subtraction, multiplication, division
- **Modular Arithmetic**: Reduction operations, cyclotomic polynomial operations
- **Coefficient Operations**: Centering, range checking, modular reduction
- **Cryptographic Support**: BigInt precision for cryptographic operations

#### Features

- `serde` - Optional serialization support
- Comprehensive test coverage
- Well-documented API with examples

#### Usage

```rust
use e3_greco_polynomial::{Polynomial, BigInt};

// Create polynomials
let poly1 = Polynomial::new(vec![BigInt::from(1), BigInt::from(2), BigInt::from(3)]);
let poly2 = Polynomial::new(vec![BigInt::from(1), BigInt::from(1)]);

// Perform operations
let sum = poly1.add(&poly2);
let product = poly1.mul(&poly2);

// Modular operations
let modulus = BigInt::from(7);
let reduced = poly1.reduce_and_center(&modulus);
```

### 🔧 e3-greco-generator

A library and CLI tool for generating BFV parameters and validation vectors:

#### Library Structure

- **`validation`**: Input validation vectors computation
- **`bounds`**: Cryptographic parameter bounds checking
- **`serialization`**: TOML/JSON serialization for Noir circuits
- **`parameters`**: BFV parameter generation

#### CLI Usage

```bash
# Generate parameters and validation vectors
cargo run --bin generator

# Output files:
# - ./output/Prover.toml (for Noir circuits)
# - ./output/validation_vectors.json
# - ./output/bounds.txt
```

#### Library Usage

```rust
use greco_generator::{
    InputValidationVectors,
    InputValidationBounds,
    ParameterGenerator,
    serialization::ProverTomlFormat,
};

// Generate BFV parameters
let params = ParameterGenerator::default_test_parameters()?;

// Compute validation vectors
let vectors = InputValidationVectors::compute(&pt, &u_rns, &e0_rns, &e1_rns, &ct, &pk)?;

// Check bounds
let bounds = InputValidationBounds::compute(&params, level)?;
bounds.check_constraints(&vectors, &zkp_modulus);

// Serialize for Noir
let toml_format = ProverTomlFormat::from_validation_vectors(&vectors);
toml_format.write_to_file("Prover.toml")?;
```

## Building and Testing

### Build All Crates

```bash
cargo build
```

### Run Tests

```bash
# Run all tests
cargo test

# Test specific crate
cargo test -p e3-greco-polynomial
cargo test -p e3-greco-generator
```

### Run CLI

```bash
# Generate parameters
cargo run --bin generator

# With release optimizations
cargo run --release --bin generator
```

## Integration with Noir Circuits

The generator produces files compatible with Noir circuits:

1. **`Prover.toml`**: Input file for Noir prover with validation vectors
2. **Constants**: Generated bounds and parameters for circuit constraints
3. **JSON**: Alternative format for integration with other tools

### Noir Integration Example

```rust
// In your Noir project's main.nr
fn main(
    pk0is: pub [Polynomial<N>; L],
    pk1is: pub [Polynomial<N>; L],
    ct0is: pub [Polynomial<N>; L],
    ct1is: pub [Polynomial<N>; L],
    u: Polynomial<N>,
    e0: Polynomial<N>,
    e1: Polynomial<N>,
    k1: Polynomial<N>,
    // ... other parameters from Prover.toml
) {
    let circuit = BfvPkEncryptionCircuit::new(/* ... */);
    circuit.correct_encryption();
}
```

## Architecture Benefits

### 🔄 Modular Design

- **Separation of Concerns**: Polynomial operations separate from parameter generation
- **Reusable Components**: Libraries can be used independently
- **Easy Testing**: Each module has focused test coverage

### 📚 Library + CLI Pattern

- **Library**: Core functionality for integration
- **CLI**: User-friendly interface for parameter generation
- **Dual Use**: Same code serves both programmatic and command-line use cases

### 🔧 Maintainable Structure

- **Clear Dependencies**: Polynomial ← Generator ← CLI
- **Documentation**: Comprehensive docs and examples
- **Error Handling**: Proper error types and propagation

## Development Workflow

### Adding New Features

1. **Polynomial Operations**: Add to `e3-greco-polynomial/src/lib.rs`
2. **Parameter Logic**: Add to appropriate module in `e3-greco-generator/src/`
3. **CLI Features**: Modify `e3-greco-generator/src/generator/main.rs`

### Testing Strategy

1. **Unit Tests**: Each function has focused tests
2. **Integration Tests**: CLI and library integration
3. **Doc Tests**: Examples in documentation are tested

### Release Process

```bash
# Test everything
cargo test

# Build release
cargo build --release

# Run CLI
./target/release/generator
```

## Future Enhancements

- [ ] Complete validation vector computation implementation
- [ ] Add more polynomial operations (GCD, factorization)
- [ ] Implement full bounds computation from BFV parameters
- [ ] Add benchmarks for performance optimization
- [ ] Support for different parameter sets
- [ ] Integration tests with actual Noir circuits

## Contributing

1. Follow the modular structure
2. Add tests for new functionality
3. Update documentation
4. Ensure CLI remains user-friendly
5. Maintain backward compatibility for library APIs
