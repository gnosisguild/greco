# Greco Generator

A Rust library and CLI tool for generating cryptographic parameters and constants for Greco zero-knowledge proofs in BFV homomorphic encryption.

## Features

- **Parameter Generation**: BFV encryption parameters and bounds computation
- **Input Validation**: Vector computation and constraint checking
- **Multiple Outputs**: Noir constants and Prover TOML file generation
- **CLI Interface**: User-friendly command-line tool with comprehensive options
- **Library API**: Can be used as a Rust crate in other projects

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
e3-greco-generator = "0.1.0"
```

## Quick Start

### CLI Usage

```bash
# Generate with default parameters
cargo run --bin generator

# Custom parameters
cargo run --bin generator -- \
  --degree 2048 \
  --plaintext-modulus 1032193 \
  --moduli "18014398492704769" \
  --output-dir my_output
```

### Library Usage

```rust
use greco_generator::{BfvConfig, GeneratorConfig, generate_all_outputs};

let bfv_config = BfvConfig {
    degree: 2048,
    plaintext_modulus: 1032193,
    moduli: vec![18014398492704769],
};

let generator_config = GeneratorConfig {
    output_dir: "output".into(),
    generate_toml: true,
};

let results = generate_all_outputs(bfv_config, generator_config)?;
```

## CLI Options

| Option                | Short | Description                               | Default             |
| --------------------- | ----- | ----------------------------------------- | ------------------- |
| `--degree`            | `-d`  | Cyclotomic polynomial degree (power of 2) | `2048`              |
| `--plaintext-modulus` | `-t`  | Plaintext modulus                         | `1032193`           |
| `--moduli`            | `-q`  | Ciphertext moduli (comma-separated)       | `18014398492704769` |
| `--output-dir`        | `-o`  | Output directory for generated files      | `output`            |
| `--no-toml`           |       | Skip generating Prover.toml file          | `false`             |

## Generated Files

### Prover TOML (`Prover.toml`)

Contains input validation vectors (plaintext coefficients, encryption randomness polynomials, ciphertext polynomials) reduced modulo the BN254 scalar field.

## Library API

### Core Types

- `BfvConfig`: BFV encryption parameters
- `GeneratorConfig`: Output generation settings
- `GenerationResults`: Generated files and computed data
- `InputValidationVectors`: Computed validation vectors
- `InputValidationBounds`: Constraint bounds

### Key Functions

- `generate_all_outputs()`: High-level generation function
- `get_default_zkp_modulus()`: Returns BN254 scalar field modulus

## License

This project is licensed under the MIT License - see the [LICENSE](../../LICENSE) file for details.
