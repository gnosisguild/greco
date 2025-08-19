# Greco Circuits Library

A Noir library implementing zero-knowledge circuits for proving correct ciphertext encryption under BFV public key homomorphic encryption.

## ⚠️ Important Disclaimer

> **Warning:** By default, this circuit enforces a binary check (only accepts `0` or `1`) for data validation.  
> If you remove that check (in `circuits/src/crypto/pk_encryption.nr`), the circuit will no longer verify that inputs are strictly `0` or `1`.  
> _Proceed with caution_: you are responsible for ensuring all inputs are valid and safe for your desired use case.

## Features

- **BFV Encryption Circuit**: Proves correct encryption of ciphertexts under BFV scheme
- **Polynomial Operations**: Efficient polynomial arithmetic in Noir
- **Safe Sponge**: Implementation of the SAFE sponge construction
- **Modular Structure**: Organized into crypto and math modules
- **Constant Parameters**: Pre-generated constants for BFV parameters

## Structure

```
circuits/
├── src/
│   ├── crypto/
│   │   ├── pk_encryption.nr   # BFV encryption circuit
│   │   └── safe.nr            # SAFE sponge implementation
│   ├── math/
│   │   └── polynomial.nr      # Polynomial arithmetic
│   └── lib.nr                 # Library entry point
```

## Usage

Add this to your `Nargo.toml`:

```toml
[dependencies]
circuits = { tag = "v0.1.0", git = "https://github.com/gnosisguild/greco" }
```

Basic usage:

```rust
use circuits::crypto::pk_encryption::BfvPkEncryptionCircuit;
use circuits::math::polynomial::Polynomial;

// Create polynomials for public keys, ciphertexts, etc.
let circuit = BfvPkEncryptionCircuit::<let N: u32, let L: u32>::new(
    pk0is, pk1is, ct0is, ct1is,
    u, e0, e1, k1,
    r1is, r2is, p1is, p2is
);

// Verify encryption correctness
circuit.correct_encryption();
```

## Constants

The constants are generated with the following parameters:

- Degree: 2048
- Plaintext modulus: 1032192
- Ciphertext moduli: [18014398492704769]

## Mathematical Background

The circuits implement zero-knowledge proofs for polynomial operations in rings of the form `Z_q[X]/(X^N + 1)`, where:

- `Z_q` is the ring of integers modulo a prime `q`
- `X^N + 1` is a cyclotomic polynomial
- Coefficients are centered in `[-(q-1)/2, (q-1)/2]`

### Circuit Components

1. **Polynomial Circuit**: Handles polynomial arithmetic and range checks
2. **SAFE Sponge**: Implements the SAFE sponge construction for generating challenge values
3. **BFV Encryption**: Proves correct encryption under the BFV scheme

## Testing

Run the test suite:

```bash
nargo test
```

## License

This project is licensed under the MIT License - see the [LICENSE](../LICENSE) file for details.
