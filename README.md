# Greco

> [!WARNING]  
> This is a research project and hasn't been audited. Use in production at your own risk.

Circuit for proving the correct encryption under the BFV fully homomorphic encryption scheme. Note that this can be generalized to any RLWE-based scheme.  
Paper: https://eprint.iacr.org/2024/594  
The implementation manages the public inputs of the circuit [in a slightly different way](https://github.com/privacy-scaling-explorations/greco/pull/30) than the one described in the paper.

## Overview

This repo includes:

- A **Noir** circuit that checks the correctness of a BFV encryption.
- A **Rust script** that generates input constants for the circuit, serializes them into TOML files, and optionally JSON.

## Requirements

- [Noirup (for installing Noir)](https://noir-lang.org/getting-started/installation/)
- [Rust toolchain](https://rustup.rs/)
- `nargo` (installed automatically via `noirup`)

### Testing Guide
The Noir circuit requires structured input values in a TOML format. These values are not hardcoded in the circuit; instead, they must be generated ahead of time and passed as input during proof generation.

To generate these inputs, we use a Rust program provided in this repository. This program mimics the BFV encryption process and outputs all the necessary values in a Noir-compatible format.

- In [rs-script](https://github.com/gnosisguild/greco/tree/noir/rs-script) run cargo with `cargo run` command to generate inputs.

- Put `pk_enc_constants_1024_2x52_2048.nr` file generated in `scripts/constants/pk_enc_constants` into the [circuits/src](https://github.com/gnosisguild/greco/tree/noir/circuits/src).

- Create a new Noir project by using `nargo new`

- Use `Prover.toml` file, which will be generated in `rs-script/scripts/pk_enc_data` as inputs to the circuit.
- Call `correct_encryption` function in main. Call public and private inputs directly to the main function similar to this:
 `fn main(pub pk0is, pub pk1is, pub ct0is, pub ct1is, u, e0, e1, k1, r1is, r2is, p1is, p2is)`

 **Follow these steps in the terminal to generate a proof with Barretenberg:**
 1. `nargo execute` 
 Compiles and executes the Noir program. Execution will generate *`./target/witness_name.gz`* and compilation of Noir program will generate *`./target/new_project_name.json`*.
 2. `bb prove -b ./target/new_project_name.json -w ./target/new_project_name.gz -o ./target`
 Generate a proof by using `Barretenberg` backend to `./target`.
 
 3. `bb write_vk -b ./target/new_project_name.json -o ./target` 
 Generate the verification key and save to `./target/vk`
 4. `bb verify -k ./target/vk -p ./target/proof`
 Verify the proof.

To create a solidity contract follow this command at the `3.` step:

3. **`bb write_vk -b ./target/new_project_name.json -o ./target --oracle_hash keccak`** --- Generate the verification key. Need to pass the `--oracle_hash keccak` flag when generating vkey and proving to instruct bb to use keccak as the hash function which is more optimal in Solidity

4. **`bb write_solidity_verifier -k ./target/vk -o ./target/Verifier.sol`**
 Generate the Solidity verifier from the vkey

For further steps to deploy and test the contract please visit [Noir Document - How to Solidity Verifier](https://noir-lang.org/docs/how_to/how-to-solidity-verifier)

