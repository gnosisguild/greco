//! Prover TOML generator.
//!
//! This module generates `Prover.toml` files containing input validation vectors
//! and BFV parameter bounds for use with Noir provers.

use crate::bounds::InputValidationBounds;
use crate::utils::to_string_1d_vec;
use crate::vectors::InputValidationVectors;
use serde::Serialize;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Generator for Prover TOML files
pub struct TomlGenerator;

/// Parameter bounds to include in the TOML
#[derive(Serialize)]
struct ProverParamsTable {
    e_bound: String,
    u_bound: String,
    k0is: Vec<String>,
    k1_low_bound: String,
    k1_up_bound: String,
    p1_bounds: Vec<String>,
    p2_bounds: Vec<String>,
    pk_bounds: Vec<String>,
    q_mod_t: String,
    qis: Vec<String>,
    r1_low_bounds: Vec<String>,
    r1_up_bounds: Vec<String>,
    r2_bounds: Vec<String>,
    tag: String,
}

/// Structure for individual vector tables in TOML
#[derive(Serialize)]
struct ProverVectorsTable {
    coefficients: Vec<String>,
}

/// Complete `Prover.toml` format including params and vectors
#[derive(Serialize)]
struct ProverTomlFormat {
    params: ProverParamsTable,
    // Original structured arrays (for backward compatibility)
    ct0is: Vec<ProverVectorsTable>,
    ct1is: Vec<ProverVectorsTable>,
    pk0is: Vec<ProverVectorsTable>,
    pk1is: Vec<ProverVectorsTable>,
    r1is: Vec<ProverVectorsTable>,
    r2is: Vec<ProverVectorsTable>,
    p1is: Vec<ProverVectorsTable>,
    p2is: Vec<ProverVectorsTable>,
    u: ProverVectorsTable,
    e0: ProverVectorsTable,
    e1: ProverVectorsTable,
    k1: ProverVectorsTable,
    // Flattened arrays for optimized circuit
    flattened_pk0is: ProverVectorsTable,
    flattened_pk1is: ProverVectorsTable,
    flattened_ct0is: ProverVectorsTable,
    flattened_ct1is: ProverVectorsTable,
    flattened_r1is: ProverVectorsTable,
    flattened_r2is: ProverVectorsTable,
    flattened_p1is: ProverVectorsTable,
    flattened_p2is: ProverVectorsTable,
}

impl TomlGenerator {
    /// Create a new TOML generator
    pub fn new() -> Self {
        Self
    }

    /// Generate `Prover.toml` file with bounds and vectors
    pub fn generate(
        &self,
        bounds: &InputValidationBounds,
        vectors: &InputValidationVectors,
        output_dir: &Path,
    ) -> Result<PathBuf, Box<dyn std::error::Error>> {
        let output_path = output_dir.join("Prover.toml");
        let mut file = File::create(&output_path)?;

        // Aggregate parameters and vectors
        let toml_data = self.to_prover_toml_format(bounds, vectors);

        // Serialize to TOML
        let toml_string = toml::to_string(&toml_data)?;

        // Write to file
        file.write_all(toml_string.as_bytes())?;

        Ok(output_path)
    }

    /// Convert bounds and vectors to `ProverTomlFormat`
    fn to_prover_toml_format(
        &self,
        bounds: &InputValidationBounds,
        vecs: &InputValidationVectors,
    ) -> ProverTomlFormat {
        ProverTomlFormat {
            params: ProverParamsTable {
                e_bound: bounds.e_bound.to_string(),
                u_bound: bounds.u_bound.to_string(),
                k0is: bounds.k0is.iter().map(|b| b.to_string()).collect(),
                k1_low_bound: bounds.k1_low_bound.to_string(),
                k1_up_bound: bounds.k1_up_bound.to_string(),
                p1_bounds: bounds.p1_bounds.iter().map(|b| b.to_string()).collect(),
                p2_bounds: bounds.p2_bounds.iter().map(|b| b.to_string()).collect(),
                pk_bounds: bounds.pk_bounds.iter().map(|b| b.to_string()).collect(),
                q_mod_t: bounds.q_mod_t.to_string(),
                qis: bounds.moduli.iter().map(|b| b.to_string()).collect(),
                r1_low_bounds: bounds.r1_low_bounds.iter().map(|b| b.to_string()).collect(),
                r1_up_bounds: bounds.r1_up_bounds.iter().map(|b| b.to_string()).collect(),
                r2_bounds: bounds.r2_bounds.iter().map(|b| b.to_string()).collect(),
                tag: bounds.tag.to_string(),
            },
            ct0is: vecs
                .ct0is
                .iter()
                .map(|v| ProverVectorsTable {
                    coefficients: to_string_1d_vec(v),
                })
                .collect(),
            ct1is: vecs
                .ct1is
                .iter()
                .map(|v| ProverVectorsTable {
                    coefficients: to_string_1d_vec(v),
                })
                .collect(),
            pk0is: vecs
                .pk0is
                .iter()
                .map(|v| ProverVectorsTable {
                    coefficients: to_string_1d_vec(v),
                })
                .collect(),
            pk1is: vecs
                .pk1is
                .iter()
                .map(|v| ProverVectorsTable {
                    coefficients: to_string_1d_vec(v),
                })
                .collect(),
            r1is: vecs
                .r1is
                .iter()
                .map(|v| ProverVectorsTable {
                    coefficients: to_string_1d_vec(v),
                })
                .collect(),
            r2is: vecs
                .r2is
                .iter()
                .map(|v| ProverVectorsTable {
                    coefficients: to_string_1d_vec(v),
                })
                .collect(),
            p1is: vecs
                .p1is
                .iter()
                .map(|v| ProverVectorsTable {
                    coefficients: to_string_1d_vec(v),
                })
                .collect(),
            p2is: vecs
                .p2is
                .iter()
                .map(|v| ProverVectorsTable {
                    coefficients: to_string_1d_vec(v),
                })
                .collect(),
            u: ProverVectorsTable {
                coefficients: to_string_1d_vec(&vecs.u),
            },
            e0: ProverVectorsTable {
                coefficients: to_string_1d_vec(&vecs.e0),
            },
            e1: ProverVectorsTable {
                coefficients: to_string_1d_vec(&vecs.e1),
            },
            k1: ProverVectorsTable {
                coefficients: to_string_1d_vec(&vecs.k1),
            },
            flattened_pk0is: ProverVectorsTable {
                coefficients: to_string_1d_vec(&vecs.flattened_pk0is()),
            },
            flattened_pk1is: ProverVectorsTable {
                coefficients: to_string_1d_vec(&vecs.flattened_pk1is()),
            },
            flattened_ct0is: ProverVectorsTable {
                coefficients: to_string_1d_vec(&vecs.flattened_ct0is()),
            },
            flattened_ct1is: ProverVectorsTable {
                coefficients: to_string_1d_vec(&vecs.flattened_ct1is()),
            },
            flattened_r1is: ProverVectorsTable {
                coefficients: to_string_1d_vec(&vecs.flattened_r1is()),
            },
            flattened_r2is: ProverVectorsTable {
                coefficients: to_string_1d_vec(&vecs.flattened_r2is()),
            },
            flattened_p1is: ProverVectorsTable {
                coefficients: to_string_1d_vec(&vecs.flattened_p1is()),
            },
            flattened_p2is: ProverVectorsTable {
                coefficients: to_string_1d_vec(&vecs.flattened_p2is()),
            },
        }
    }
}
