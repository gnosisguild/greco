use circuits::crypto::{Params, pk_encryption::BfvPkEncryptionCircuit};
use circuits::math::polynomial::Polynomial;

fn main(
    params: pub Params<2048, 1>,
    pk0is: pub [Polynomial<2048>; 1],
    pk1is: pub [Polynomial<2048>; 1],
    ct0is: pub [Polynomial<2048>; 1],
    ct1is: pub [Polynomial<2048>; 1],
    u: Polynomial<2048>,
    e0: Polynomial<2048>,
    e1: Polynomial<2048>,
    k1: Polynomial<2048>,
    r1is: [Polynomial<4095>; 1],
    r2is: [Polynomial<2047>; 1],
    p1is: [Polynomial<4095>; 1],
    p2is: [Polynomial<2047>; 1],
) {
    let circuit = BfvPkEncryptionCircuit::new(
        params,
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
        p2is,
    );
    circuit.correct_encryption();
}
