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

- [Noirup (for installing Noir)](https://noir-lang.org/docs/getting_started/quick_start#noir)
- [BBup (for installing Barretenberg)](https://noir-lang.org/docs/getting_started/quick_start#proving-backend)
- [Rust toolchain](https://rustup.rs/)
- `nargo` (installed automatically via `noirup`)

---

## Testing Guide

The Noir circuit requires structured input values in a TOML format. These values are not hardcoded in the circuit; instead, they must be generated ahead of time and passed as input during proof generation.

To generate these inputs, we use a Rust program provided in this repository. This program mimics the BFV encryption process and outputs all necessary values in a Noir-compatible format.

- Copy Greco files into a new file. Lets name it as Noir-Greco
- In the [rs-script](https://github.com/gnosisguild/greco/tree/noir/rs-script) folder, run `cargo run` to generate inputs and constants.
- Copy the `pk_enc_constants_1024_2x52_2048.nr` file generated in `greco/rs-script/scripts/constants/pk_enc_constants` into the [greco/circuits/src](https://github.com/gnosisguild/greco/tree/noir/circuits/src) directory.
- Create a new Noir project by using command below in Noir-Greco folder
   ```rust
   nargo new <new_project_name>
   ```
   This command will create two files,

   `Noir-Greco/<new_project_name>/src/main.nr`

   `Noir-Greco/<new_project_name>/Nargo.toml`
- Open Nargo.toml and add this path to the dependency
   ```rust
   [dependencies]
   circuits = { path = "../greco/circuits", type = "lib" }
   ```
- In `main.nr`, call the `correct_encryption` function and pass public and private inputs directly in the function signature, for example:

   ```rust
   fn main(
      pk0is: pub [Polynomial<N>; L],
      pk1is: pub [Polynomial<N>; L],
      ct0is: pub [Polynomial<N>; L],
      ct1is: pub [Polynomial<N>; L],
      u: Polynomial<N>,
      e0: Polynomial<N>,
      e1: Polynomial<N>,
      k1: Polynomial<N>,
      r1is: [Polynomial<(2 * N) - 1>; L],
      r2is: [Polynomial<N - 1>; L],
      p1is: [Polynomial<(2 * N) - 1>; L],
      p2is: [Polynomial<N - 1>; L]
   )  {
      let circuit = BfvPkEncryptionCircuit::new(
         pk0is,
         pk1is,
         ct0is,
         ct1is,
         u,
         e0,
         e1,
         k1,
         r1is,
         r2is,
         p1is,
         p2is
      );
      circuit.correct_encryption();
   }
   ```

- Move `Prover.toml` file (generated in `greco/rs-script/scripts/pk_enc_data` as input to the circuit) to the `Noir-Greco/<new_project_name>`.

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
