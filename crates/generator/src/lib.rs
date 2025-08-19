//! # Greco Generator Library
//!
//! A library for generating cryptographic parameters and constants for Greco zero-knowledge proofs
//! in the context of BFV homomorphic encryption.

pub mod bfv;
pub mod bounds;
pub mod cli;
pub mod generators;
pub mod utils;
pub mod vectors;

// Re-export main types that currently work
pub use bfv::{BfvConfig, BfvHelper, EncryptionData};
pub use bounds::InputValidationBounds;
pub use cli::CliConfig;
pub use generators::toml::TomlGenerator;
pub use vectors::InputValidationVectors;

use num_traits::Num;
use polynomial::BigInt;
use std::path::PathBuf;

/// Configuration for output generation
#[derive(Clone, Debug)]
pub struct GeneratorConfig {
    pub output_dir: PathBuf,
    pub generate_toml: bool,
}

impl Default for GeneratorConfig {
    fn default() -> Self {
        Self {
            output_dir: "output".into(),
            generate_toml: true,
        }
    }
}

/// Results from generation process
#[derive(Debug)]
pub struct GenerationResults {
    pub vectors: InputValidationVectors,
    pub bounds: InputValidationBounds,
    pub toml_file: Option<PathBuf>,
}

/// High-level function to generate all outputs given BFV configuration
pub fn generate_all_outputs(
    bfv_config: BfvConfig,
    generator_config: GeneratorConfig,
) -> Result<GenerationResults, Box<dyn std::error::Error>> {
    // Create BFV helper and generate encryption
    let helper = BfvHelper::new(bfv_config)?;
    let encryption_data = helper.generate_sample_encryption()?;

    // Compute input validation vectors
    let vectors = InputValidationVectors::compute(
        &encryption_data.plaintext,
        &encryption_data.u_rns,
        &encryption_data.e0_rns,
        &encryption_data.e1_rns,
        &encryption_data.ciphertext,
        &encryption_data.public_key,
        &helper.params,
    )?;

    // Compute bounds
    let bounds = InputValidationBounds::compute(&helper.params, encryption_data.plaintext.level())?;

    // Get ZKP modulus
    let zkp_modulus = BigInt::from_str_radix(
        "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        10,
    )?;

    // Check constraints
    bounds.check_constraints(&vectors, &zkp_modulus);

    // Create output directory
    std::fs::create_dir_all(&generator_config.output_dir)?;

    let mut results = GenerationResults {
        vectors: vectors.clone(),
        bounds: bounds.clone(),
        toml_file: None,
    };

    // Generate Prover TOML if requested
    if generator_config.generate_toml {
        let toml_generator = TomlGenerator::new();
        let toml_path = toml_generator.generate(
            &bounds,
            &vectors.standard_form(&zkp_modulus),
            &generator_config.output_dir,
        )?;
        results.toml_file = Some(toml_path);
    }

    Ok(results)
}

/// Test function to check what specific errors we get with vectors
#[cfg(test)]
#[test]
pub fn test_vectors_computation() -> Result<(), Box<dyn std::error::Error>> {
    let config = BfvConfig {
        degree: 2048,
        plaintext_modulus: 1032192,
        moduli: vec![18014398492704769],
    };

    let helper = BfvHelper::new(config)?;
    let encryption_data = helper.generate_sample_encryption()?;

    // Try to compute vectors - this will show us the exact errors
    let _vectors = InputValidationVectors::compute(
        &encryption_data.plaintext,
        &encryption_data.u_rns,
        &encryption_data.e0_rns,
        &encryption_data.e1_rns,
        &encryption_data.ciphertext,
        &encryption_data.public_key,
        &helper.params,
    )?;

    println!("Vectors computation successful!");
    Ok(())
}

/// Test bounds computation
#[cfg(test)]
#[test]
pub fn test_bounds_computation() -> Result<(), Box<dyn std::error::Error>> {
    let config = BfvConfig {
        degree: 2048,
        plaintext_modulus: 1032192,
        moduli: vec![18014398492704769],
    };

    let helper = BfvHelper::new(config)?;
    let encryption_data = helper.generate_sample_encryption()?;

    // Try to compute bounds
    let _bounds =
        InputValidationBounds::compute(&helper.params, encryption_data.plaintext.level())?;

    println!("Bounds computation successful!");
    Ok(())
}
