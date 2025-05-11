# Greco

> [!WARNING]  
> This is a research project and hasn't been audited. Use in production at your own risk.

Circuit for proving the correct encryption under the BFV fully homomorphic encryption scheme. Note that this can be generalized to any RLWE-based scheme.  
Paper: https://eprint.iacr.org/2024/594  
The implementation manages the public inputs of the circuit [in a slightly different way](https://github.com/privacy-scaling-explorations/greco/pull/30) than the one described in the paper.

## Overview

This repo includes:

- A **Noir** circuit that checks the correctness of BFV encryption.
- A **Rust script** that generates input constants for the circuit, serialized into TOML files (and optionally JSON).

## Requirements

- [Noirup (for installing Noir)](https://noir-lang.org/getting-started/installation/)
- [Rust toolchain](https://rustup.rs/)
- `nargo` (installed automatically via `noirup`)

---

## Testing Guide

The Noir circuit requires structured input values in a TOML format. These values are not hardcoded in the circuit; instead, they must be generated ahead of time and passed as input during proof generation.

To generate these inputs, we use a Rust program provided in this repository. This program mimics the BFV encryption process and outputs all necessary values in a Noir-compatible format.

- In the [rs-script](https://github.com/gnosisguild/greco/tree/noir/rs-script) folder, run `cargo run` to generate inputs and constants.
- Copy the `pk_enc_constants_1024_2x52_2048.nr` file generated in `scripts/constants/pk_enc_constants` into the [circuits/src](https://github.com/gnosisguild/greco/tree/noir/circuits/src) directory.
- Create a new Noir project using `nargo new`.
- Use the `Prover.toml` file (generated in `rs-script/scripts/pk_enc_data`) as input to the circuit.
- In `main.nr`, call the `correct_encryption` function and pass public and private inputs directly in the function signature, for example:

```rust
fn main(pub pk0is, pub pk1is, pub ct0is, pub ct1is, u, e0, e1, k1, r1is, r2is, p1is, p2is)
```

---

### Generating a Proof with Barretenberg

1. **Compile and execute the Noir program**  
   ```bash
   nargo execute
   ```
   This command compiles and executes the circuit, generating `./target/<witness_name>.gz` and `./target/<new_project_name>.json`.

2. **Generate a proof**  
   ```bash
   bb prove -b ./target/new_project_name.json -w ./target/new_project_name.gz -o ./target
   ```
   This command generates a proof using the Barretenberg backend.

3. **Write the verification key**  
   ```bash
   bb write_vk -b ./target/new_project_name.json -o ./target
   ```
   This command generates the verification key.

4. **Verify the proof**  
   ```bash
   bb verify -k ./target/vk -p ./target/proof
   ```
   This command verifies the generated proof.

---

### Creating a Solidity Verifier

After completing step 2 above (generating the proof), you can proceed with:

3. **Generate the verification key using Keccak**  
   ```bash
   bb write_vk -b ./target/new_project_name.json -o ./target --oracle_hash keccak
   ```
   This command generates the verification key using Keccak. You must pass the `--oracle_hash keccak` flag when generating the vkey and during proving to ensure compatibility with Solidity-based verification.

4. **Generate the Solidity verifier**  
   ```bash
   bb write_solidity_verifier -k ./target/vk -o ./target/Verifier.sol
   ```
   This command generates a Solidity verifier from the verification key.

For further steps to deploy and test the contract, refer to the [Noir documentation on Solidity verifiers](https://noir-lang.org/docs/how_to/how-to-solidity-verifier).
