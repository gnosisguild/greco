use std::fs::File;
use std::io::Write;

use axiom_eth::rlc::{
    circuit::{builder::RlcCircuitBuilder, instructions::RlcCircuitInstructions},
    utils::executor::RlcExecutor,
};

use greco::pk_encryption_circuit::BfvPkEncryptionCircuit;
use halo2_base::{
    gates::circuit::CircuitBuilderStage,
    halo2_proofs::{
        halo2curves::bn256::{Bn256, Fr, G1Affine},
        plonk::{keygen_pk, keygen_vk, ProvingKey, VerifyingKey},
        poly::kzg::commitment::ParamsKZG,
        SerdeFormat,
    },
    utils::fs::gen_srs,
};
use halo2_solidity_verifier::{BatchOpenScheme::Bdfg21, SolidityGenerator};
use serde::{Deserialize, Serialize};

// ----------------------
// Key Storage Structure
// ----------------------

#[derive(Serialize, Deserialize)]
struct EncryptionKeys {
    proving_key: Vec<u8>,
}

// ----------------------
// Main Execution
// ----------------------

fn main() {
    // Initialize circuit parameters
    let empty_circuit = BfvPkEncryptionCircuit::create_empty_circuit(1, 2048);
    let public_instances: Vec<Vec<Fr>> = empty_circuit.instances();

    // Generate structured reference string
    let srs_degree = 17; // k
    let kzg_params = gen_srs(srs_degree);

    // Configure and build RLC circuit
    let mut circuit_builder = RlcCircuitBuilder::<Fr>::from_stage(CircuitBuilderStage::Keygen, 0)
        .use_k(srs_degree as usize);

    circuit_builder.base.set_lookup_bits((srs_degree - 1) as usize);
    circuit_builder.base.set_instance_columns(1);

    let rlc_executor = RlcExecutor::new(circuit_builder, empty_circuit.clone());
    rlc_executor.0.calculate_params(Some(9));

    // Generate cryptographic keys
    let (verification_key, proving_key) = generate_zk_keys(&kzg_params, &rlc_executor);

    // Generate verifier contract
    let solidity_verifier =
        create_solidity_verifier(&kzg_params, &verification_key, public_instances[0].len());

    // Persist generated artifacts
    save_keys_to_file(&proving_key).expect("Failed to save keys to file");

    save_solidity_verifier(&solidity_verifier).expect("Failed to save Solidity verifier");

    println!(
        "Successfully generated and stored:\n\
        - Cryptographic keys: bfv_encryption_keys.json\n\
        - Verifier contract: BfvEncryptionVerifier.sol"
    );
}

// ----------------------
// Key Generation Logic
// ----------------------

fn generate_zk_keys(
    kzg_params: &ParamsKZG<Bn256>,
    circuit: &axiom_eth::rlc::utils::two_phase::TwoPhaseCircuit<
        Fr,
        RlcExecutor<Fr, BfvPkEncryptionCircuit>,
    >,
) -> (VerifyingKey<G1Affine>, ProvingKey<G1Affine>) {
    let verification_key: axiom_eth::halo2_proofs::plonk::VerifyingKey<
        axiom_eth::halo2curves::bn256::G1Affine,
    > = keygen_vk(kzg_params, circuit).expect("Failed to generate verification key");

    let proving_key = keygen_pk(kzg_params, verification_key.clone(), circuit)
        .expect("Failed to generate proving key");

    (verification_key, proving_key)
}

// ----------------------
// Verifier Generation
// ----------------------

fn create_solidity_verifier(
    kzg_params: &ParamsKZG<Bn256>,
    verification_key: &VerifyingKey<G1Affine>,
    num_inputs: usize,
) -> String {
    SolidityGenerator::new(kzg_params, verification_key, Bdfg21, num_inputs)
        .render()
        .expect("Failed to generate Solidity verifier")
}

// ----------------------
// File I/O Operations
// ----------------------

fn save_keys_to_file(
    proving_key: &ProvingKey<G1Affine>,
) -> std::io::Result<()> {
    let keys = EncryptionKeys {
        proving_key: proving_key.to_bytes(SerdeFormat::Processed),
    };

    let keys_json = serde_json::to_string_pretty(&keys).expect("Failed to serialize keys to JSON");

    File::create("./keys/bfv_encryption_keys.json")?.write_all(keys_json.as_bytes())
}

fn save_solidity_verifier(contract_code: &str) -> std::io::Result<()> {
    File::create("./contracts/BfvEncryptionVerifier.sol")?.write_all(contract_code.as_bytes())
}
