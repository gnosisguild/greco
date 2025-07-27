use circuits::crypto::{Params, pk_encryption::BfvPkEncryptionCircuit};
use circuits::math::polynomial::Polynomial;

fn main(
    params: pub Params<8192, 3>,
    pk0is: pub [Polynomial<8192>; 3],
    pk1is: pub [Polynomial<8192>; 3],
    ct0is: pub [Polynomial<8192>; 3],
    ct1is: pub [Polynomial<8192>; 3],
    u: Polynomial<8192>,
    e0: Polynomial<8192>,
    e1: Polynomial<8192>,
    k1: Polynomial<8192>,
    r1is: [Polynomial<16383>; 3],
    r2is: [Polynomial<8191>; 3],
    p1is: [Polynomial<16383>; 3],
    p2is: [Polynomial<8191>; 3],
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
