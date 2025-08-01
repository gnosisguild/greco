use circuits::crypto::{Params, pk_encryption::BfvPkEncryptionCircuit};
use circuits::math::polynomial::Polynomial;

fn main(
    params: pub Params<128, 2>,
    pk0is: pub [Polynomial<128>; 2],
    pk1is: pub [Polynomial<128>; 2],
    ct0is: pub [Polynomial<128>; 2],
    ct1is: pub [Polynomial<128>; 2],
    u: Polynomial<128>,
    e0: Polynomial<128>,
    e1: Polynomial<128>,
    k1: Polynomial<128>,
    r1is: [Polynomial<255>; 2],
    r2is: [Polynomial<127>; 2],
    p1is: [Polynomial<255>; 2],
    p2is: [Polynomial<127>; 2],
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
