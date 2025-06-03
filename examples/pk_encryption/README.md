# Public Key Encryption Example

This example demonstrates a zero-knowledge proof for BFV public key encryption operations using Noir. The circuit proves the correctness of ciphertext generation using a public key while keeping the plaintext and encryption randomness private.

## Dependencies

The `main.nr` circuit implements the following using core `circuits` library.

- Verification of BFV ciphertext components
- Polynomial arithmetic with proper modular reduction
- Bounds checking for encryption randomness

The circuit takes several inputs:

- `pk0is`, `pk1is`: Public key components
- `ct0is`, `ct1is`: Ciphertext components
- `u`, `e0`, `e1`: Encryption randomness
- `k1`: Additional parameter
- `r1is`, `r2is`, `p1is`, `p2is`: Polynomial reduction parameters

## Running the Example

1. Make sure you have Noir and Barretenberg installed (see main [README](../README.md) for installation instructions)

2. From the `pk_encryption` directory, compile the circuit:

```bash
nargo compile
```

3. To generate valid inputs for the circuit:

From the project root, run:

```bash
cargo run --bin generator -- -o ../examples/pk_encryption/Prover.toml && rm -rf ../examples/pk_encryption/constants.nr
```

This will create a `Prover.toml` file with the necessary polynomial inputs formatted for the circuit (we remove the `constants.nr` since they are not necessary).

4. Execute the circuit:

```bash
nargo execute
```

5. Generate the proof:

```bash
bb prove -b ./target/pk_encryption.json -w ./target/pk_encryption.gz -o ./target
```

5. Verify the proof:

```bash
bb verify -k ./target/vk -p ./target/proof
```
