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
- Noir circuit constants generation
- Prover TOML file generation

```rust
use e3_greco_generator::{BfvConfig, GeneratorConfig, generate_all_outputs};

let config = BfvConfig {
    degree: 1024,
    plaintext_modulus: 2048,
    moduli: vec![4503599625535489, 4503599626321921],
};

let results = generate_all_outputs(config, GeneratorConfig::default())?;
```

## Quick Start

### Prerequisites

1. Install Rust:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

2. Install Noir using noirup:

```bash
curl -L https://raw.githubusercontent.com/noir-lang/noirup/refs/heads/main/install | bash
noirup
```

3. Install Barretenberg (BB) proving backend:

```bash
curl -L https://raw.githubusercontent.com/AztecProtocol/aztec-packages/refs/heads/master/barretenberg/bbup/install | bash
bbup
```

### Building and Testing

1. Build and test Rust components:

```bash
# Build all Rust components
cargo build

# Run Rust tests
cargo test

# Check Rust formatting
cargo fmt --check
```

2. Build and test Noir circuits:

```bash
# Check Noir formatting
nargo fmt --check

# Build Noir circuits
nargo check

# Run Noir tests
nargo test
```

3. Generate parameters:

```bash
cargo run --bin generator
```

This will create:

- `constants.nr`: Noir circuit parameters and bounds
- `Prover.toml`: Input validation vectors

## Mathematical Background

The workspace implements polynomial arithmetic in rings of the form `Z_q[X]/(X^N + 1)`, where:

- `Z_q` is the ring of integers modulo a prime `q`
- `X^N + 1` is a cyclotomic polynomial
- Coefficients are centered in `[-(q-1)/2, (q-1)/2]`

These structures are fundamental to the BFV homomorphic encryption scheme and its zero-knowledge proofs.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
