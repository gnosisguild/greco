use circuits::crypto::{Params, pk_encryption::BfvPkEncryptionCircuit};
use circuits::math::polynomial::Polynomial;

fn main(
    params: pub Params<1024, 1>,
    // Use the flattened arrays that are actually in the TOML
    flattened_pk0is: pub Polynomial<1024>,
    flattened_pk1is: pub Polynomial<1024>,
    flattened_ct0is: pub Polynomial<1024>,
    flattened_ct1is: pub Polynomial<1024>,
    u: Polynomial<1024>,
    e0: Polynomial<1024>,
    e1: Polynomial<1024>,
    k1: Polynomial<1024>,
    flattened_r1is: Polynomial<2047>,
    flattened_r2is: Polynomial<1023>,
    flattened_p1is: Polynomial<2047>,
    flattened_p2is: Polynomial<1023>,
) {
    let circuit = BfvPkEncryptionCircuit::new(
        params,
        flattened_pk0is.coefficients,
        flattened_pk1is.coefficients,
        flattened_ct0is.coefficients,
        flattened_ct1is.coefficients,
        u,
        e0,
        e1,
        k1,
        flattened_r1is.coefficients,
        flattened_r2is.coefficients,
        flattened_p1is.coefficients,
        flattened_p2is.coefficients,
    );
    circuit.correct_encryption();
}
