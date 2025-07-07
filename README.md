# Greco

A Rust workspace for generating and validating zero-knowledge proofs of correct ciphertext encryption under BFV public key homomorphic encryption.

> This project is based on the research and implementation by the [Privacy & Scaling Explorations (PSE) team](https://github.com/privacy-scaling-explorations/greco). We extend our gratitude for their groundbreaking work on zero-knowledge proofs for BFV encryption correctness, detailed in their [research paper](https://eprint.iacr.org/2024/594).

## Crates

### `e3-greco-polynomial`

Core polynomial arithmetic library supporting:

- Basic operations (add, sub, mul, div) with arbitrary precision
- Modular reduction and centered coefficients
- Cyclotomic polynomial operations
- Range checking for cryptographic bounds

```rust
use e3_greco_polynomial::{Polynomial, BigInt};

let poly = Polynomial::new(vec![BigInt::from(2), BigInt::from(3), BigInt::from(1)]); // 2x² + 3x + 1
let reduced = poly.reduce_and_center(&modulus); // Center coefficients in [-q/2, q/2]
```

### `e3-greco-generator`

Generator for cryptographic parameters and constants:

- BFV parameter generation and validation
- Input validation vector computation
- Prover TOML file generation

```rust
use e3_greco_generator::{BfvConfig, GeneratorConfig, generate_all_outputs};

let config = BfvConfig {
    degree: 2048,
    plaintext_modulus: 1032193,
    moduli: vec![4503599626321921],
};

let results = generate_all_outputs(config, GeneratorConfig::default())?;
```

## Quick Start

1. Install dependencies:

```bash
cargo build --workspace
```

2. Generate parameters:

```bash
cargo run --bin generator
```

This will create:

- `Prover.toml`: Input validation vectors

## Mathematical Background

The workspace implements polynomial arithmetic in rings of the form `Z_q[X]/(X^N + 1)`, where:

- `Z_q` is the ring of integers modulo a prime `q`
- `X^N + 1` is a cyclotomic polynomial
- Coefficients are centered in `[-(q-1)/2, (q-1)/2]`

These structures are fundamental to the BFV homomorphic encryption scheme and its zero-knowledge proofs.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
