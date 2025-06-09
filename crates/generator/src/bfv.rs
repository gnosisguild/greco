//! BFV parameter utilities and encryption helpers.
//!
//! This module provides utilities for working with BFV encryption parameters,
//! generating sample encryptions, and managing encryption contexts.

use fhe::bfv::{
    BfvParameters, BfvParametersBuilder, Ciphertext, Encoding, Plaintext, PublicKey, SecretKey,
};
use fhe_math::{rq::Poly, zq::Modulus};
use fhe_traits::*;
use num_bigint::BigInt;
use num_traits::Num;
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::sync::Arc;

/// Configuration for BFV parameters
#[derive(Clone, Debug)]
pub struct BfvConfig {
    pub degree: usize,
    pub plaintext_modulus: u64,
    pub moduli: Vec<u64>,
}

impl Default for BfvConfig {
    fn default() -> Self {
        Self {
            degree: 2048,
            plaintext_modulus: 1032193,
            moduli: vec![18014398492704769],
        }
    }
}

/// Data from a sample BFV encryption
pub struct EncryptionData {
    pub plaintext: Plaintext,
    pub ciphertext: Ciphertext,
    pub public_key: PublicKey,
    pub secret_key: SecretKey,
    pub u_rns: Poly,
    pub e0_rns: Poly,
    pub e1_rns: Poly,
}

/// Helper for working with BFV parameters and operations
pub struct BfvHelper {
    pub params: Arc<BfvParameters>,
}

impl BfvHelper {
    /// Create a new BFV helper from configuration
    pub fn new(config: BfvConfig) -> Result<Self, Box<dyn std::error::Error>> {
        let params = BfvParametersBuilder::new()
            .set_degree(config.degree)
            .set_plaintext_modulus(config.plaintext_modulus)
            .set_moduli(&config.moduli)
            .build_arc()?;

        Ok(BfvHelper { params })
    }

    /// Generate a sample encryption with all the data needed for input validation
    pub fn generate_sample_encryption(&self) -> Result<EncryptionData, Box<dyn std::error::Error>> {
        let mut rng = StdRng::seed_from_u64(0);

        // Generate keys
        let sk = SecretKey::random(&self.params, &mut rng);
        let pk = PublicKey::new(&sk, &mut rng);

        // Create a sample plaintext with some random values, in here we are assiging 3 to all the
        // coefficients
        let mut message_data = vec![3u64; self.params.degree() as usize];

        //For Crisp, the user casts the vote in the right coefficient (message_data[0]). A vote is
        //a value in {0,1}. Any other value will result in a proof that will be rejected by the Verifier.
        message_data[0] = 1;

        let pt = Plaintext::try_encode(&message_data, Encoding::poly(), &self.params)?;

        // Use extended encryption to get the polynomial data
        let (_ct, u_rns, e0_rns, e1_rns) = pk.try_encrypt_extended(&pt, &mut rng)?;

        Ok(EncryptionData {
            plaintext: pt,
            ciphertext: _ct,
            public_key: pk,
            secret_key: sk,
            u_rns,
            e0_rns,
            e1_rns,
        })
    }

    /// Get the default ZKP modulus
    pub fn default_zkp_modulus() -> BigInt {
        BigInt::from_str_radix(
            "21888242871839275222246405745257275088548364400416034343698204186575808495617",
            10,
        )
        .expect("Invalid ZKP modulus")
    }
}
