
use greco_script::*;
use fhe::bfv::{
     BfvParametersBuilder, Encoding, Plaintext, PublicKey, SecretKey,
};
use fhe_math::zq::Modulus;
use fhe_traits::*;

use num_bigint::BigInt;
use num_traits::Num;
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::path::Path;
use std::vec;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Set up the BFV parameters

    let N: u64 = 1024;

    // let plaintext_modulus: u64 = 65537;
    let plaintext_modulus: u64 = 2048;

    // let moduli: Vec<u64> = vec![
    //   1152921504606584833,
    //   1152921504598720513,
    //   1152921504597016577,
    //   1152921504595968001,
    //   1152921504595640321,
    //   1152921504593412097,
    //   1152921504592822273,
    //   1152921504592429057,
    //   1152921504589938689,
    //   1152921504586530817,
    //   1152921504585547777,
    //   1152921504583647233,
    //   1152921504581877761,
    //   1152921504581419009,
    //   1152921504580894721
    // ];

    let moduli: Vec<u64> = vec![4503599625535489, 4503599626321921];
    //let moduli: Vec<u64> = vec![4503599625535489, 4503599626321921];
    //let moduli: Vec<u64> = vec![1038337,18014398492704769,4503599625535489, 4503599626321921];
    //let moduli: Vec<u64> = vec![1038337];
    //let moduli: Vec<u64> = vec![
    //18014398492704769
    //];

    //let moduli_sizes = [20];

    let params = BfvParametersBuilder::new()
        .set_degree(N as usize)
        .set_plaintext_modulus(plaintext_modulus)
        .set_moduli(&moduli)
        //.set_moduli_sizes(&moduli_sizes)
        .build_arc()?;

    // Extract plaintext modulus
    let t = Modulus::new(params.plaintext())?;

    // Use a seedable rng for experimental reproducibility
    let mut rng = StdRng::seed_from_u64(0);

    // Generate the secret and public keys
    let sk = SecretKey::random(&params, &mut rng);

    let pk = PublicKey::new(&sk, &mut rng);

    //Sample a message and encrypt
    // let m = t.random_vec(N as usize, &mut rng);
    let m = vec![1; N as usize];
    // println!("m: {:?}", m);
    //let m: Vec<i64> = (-(N as i64 / 2)..(N as i64 / 2)).collect(); // m here is from lowest degree to largest as input into fhe.rs (REQUIRED)
    let pt = Plaintext::try_encode(&m, Encoding::poly(), &params)?;

    // Extract context
    let ctx = params.ctx_at_level(pt.level())?.clone();
    let (ct, u_rns, e0_rns, e1_rns) = pk.try_encrypt_extended(&pt, &mut rng)?;

    // Sanity check. m = Decrypt(ct)

    //let m_decrypted = unsafe { t.center_vec_vt(&sk.try_decrypt(&ct)?.value.into_vec()) };
    let m_decrypted = sk.try_decrypt(&ct)?.value.into_vec();
    assert_eq!(m_decrypted, m);

    let moduli = ctx.moduli();
    // Initialize zk proving modulus
    let p = BigInt::from_str_radix(
        "21888242871839275222246405745257275088548364400416034343698204186575808495617",
        10,
    )?;

    // Compute input validation vectors
    let res = InputValidationVectors::compute(&pt, &u_rns, &e0_rns, &e1_rns, &ct, &pk)?;
    // return Ok(());

    // Create output json with standard form polynomials
    let json_data = res.standard_form(&p).to_json();

    // Calculate bounds ---------------------------------------------------------------------
    let bounds = InputValidationBounds::compute(&params, pt.level())?;

    // Check the constraints
    bounds.check_constraints(&res, &p);

    let moduli_bitsize = {
        if let Some(&max_value) = ctx.moduli().iter().max() {
            64 - max_value.leading_zeros()
        } else {
            0
        }
    };

    // Write out files ----------------------------------------------------------------------
    let output_path = Path::new("src").join("data").join("pk_enc_data");

    // Generate filename and write file
    let filename = format!(
        "pk_enc_{}_{}x{}_{}.json",
        N,
        moduli.len(),
        moduli_bitsize,
        t.modulus()
    );
    write_json_to_file(&output_path, &filename, &json_data);

    // Generate zeros filename and write file
    let filename_zeroes = format!(
        "pk_enc_{}_{}x{}_{}_zeroes.json",
        N,
        moduli.len(),
        moduli_bitsize,
        t.modulus()
    );
    let zeroes_json = InputValidationVectors::new(moduli.len(), params.degree()).to_json();
    write_json_to_file(&output_path, &filename_zeroes, &zeroes_json);

    let filename_constants = format!(
        "pk_enc_constants_{}_{}x{}_{}.rs",
        N,
        moduli.len(),
        moduli_bitsize,
        t.modulus()
    );
    bounds.to_file(&params, &filename_constants)?;

    Ok(())
}
