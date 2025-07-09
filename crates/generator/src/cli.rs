//! CLI interface for the Greco generator.
//!
//! This module provides command-line argument parsing and the main CLI entry point
//! for the generator binary.

use crate::{generate_all_outputs, BfvConfig, GeneratorConfig};
use clap::{Arg, Command};
use std::path::PathBuf;

/// CLI configuration structure
#[derive(Debug, Clone)]
pub struct CliConfig {
    pub bfv_config: BfvConfig,
    pub generator_config: GeneratorConfig,
}

impl CliConfig {
    /// Parse CLI arguments and create configuration
    pub fn from_args() -> Result<Self, Box<dyn std::error::Error>> {
        let matches = Command::new("e3-greco-generator")
            .about(
                "Generates cryptographic parameters and constants for Greco zero-knowledge proofs",
            )
            .arg(
                Arg::new("degree")
                    .long("degree")
                    .short('d')
                    .value_name("DEGREE")
                    .help("Cyclotomic polynomial degree (power of 2)")
                    .default_value("2048"),
            )
            .arg(
                Arg::new("plaintext-modulus")
                    .long("plaintext-modulus")
                    .short('t')
                    .value_name("MODULUS")
                    .help("Plaintext modulus")
                    .default_value("1032193"),
            )
            .arg(
                Arg::new("moduli")
                    .long("moduli")
                    .short('q')
                    .value_name("MODULI")
                    .help("Ciphertext moduli (comma-separated)")
                    .default_value("18014398492704769"),
            )
            .arg(
                Arg::new("output-dir")
                    .long("output-dir")
                    .short('o')
                    .value_name("DIR")
                    .help("Output directory for generated files")
                    .default_value("output"),
            )
            .arg(
                Arg::new("no-toml")
                    .long("no-toml")
                    .help("Skip generating Prover.toml file")
                    .action(clap::ArgAction::SetTrue),
            )
            .get_matches();

        // Parse degree
        let degree: usize = matches
            .get_one::<String>("degree")
            .unwrap()
            .parse()
            .map_err(|e| format!("Invalid degree: {}", e))?;

        // Validate degree is power of 2
        if !degree.is_power_of_two() {
            return Err("Degree must be a power of 2".into());
        }

        // Parse plaintext modulus
        let plaintext_modulus: u64 = matches
            .get_one::<String>("plaintext-modulus")
            .unwrap()
            .parse()
            .map_err(|e| format!("Invalid plaintext modulus: {}", e))?;

        // Parse moduli
        let moduli_str = matches.get_one::<String>("moduli").unwrap();
        let moduli: Result<Vec<u64>, _> = moduli_str
            .split(',')
            .map(|s| s.trim().parse::<u64>())
            .collect();
        let moduli = moduli.map_err(|e| format!("Invalid moduli: {}", e))?;

        // Parse generation flags
        let generate_toml = !matches.get_flag("no-toml");

        // Parse output directory
        let output_dir = PathBuf::from(matches.get_one::<String>("output-dir").unwrap());

        let bfv_config = BfvConfig {
            degree,
            plaintext_modulus,
            moduli,
        };

        let generator_config = GeneratorConfig {
            output_dir,
            generate_toml,
        };

        Ok(CliConfig {
            bfv_config,
            generator_config,
        })
    }

    /// Run the CLI application
    pub fn run(&self) -> Result<(), Box<dyn std::error::Error>> {
        println!("Greco Generator");
        println!("Generating cryptographic parameters...");

        // Display configuration
        println!("\nConfiguration:");
        println!("- Degree: {}", self.bfv_config.degree);
        println!("- Plaintext modulus: {}", self.bfv_config.plaintext_modulus);
        println!("- Ciphertext moduli: {:?}", self.bfv_config.moduli);
        println!("- Output directory: {:?}", self.generator_config.output_dir);
        println!("- Generate TOML: {}", self.generator_config.generate_toml);

        // Generate outputs
        let results = generate_all_outputs(self.bfv_config.clone(), self.generator_config.clone())?;

        println!("\nOutputs:");
        if let Some(toml_path) = &results.toml_file {
            println!("- Prover TOML: {:?}", toml_path);
        }

        Ok(())
    }
}

/// Main CLI entry point
pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    let config = CliConfig::from_args()?;
    config.run()
}
