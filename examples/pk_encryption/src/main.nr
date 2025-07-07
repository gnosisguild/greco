use circuits::crypto::{Params, pk_encryption::BfvPkEncryptionCircuit};
use circuits::math::polynomial::Polynomial;

fn main(
    params: pub Params<1024, 2>,
    pk0is: pub [Polynomial<1024>; 2],
    pk1is: pub [Polynomial<1024>; 2],
    ct0is: pub [Polynomial<1024>; 2],
    ct1is: pub [Polynomial<1024>; 2],
    u: Polynomial<1024>,
    e0: Polynomial<1024>,
    e1: Polynomial<1024>,
    k1: Polynomial<1024>,
    r1is: [Polynomial<2047>; 2],
    r2is: [Polynomial<1023>; 2],
    p1is: [Polynomial<2047>; 2],
    p2is: [Polynomial<1023>; 2],
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
