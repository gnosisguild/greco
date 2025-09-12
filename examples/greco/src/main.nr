use circuits::crypto::greco::{Greco, Params};
use circuits::math::polynomial::Polynomial;

fn main(
    params: pub Params<8192, 4>,
    pk0is: pub [Polynomial<8192>; 4],
    pk1is: pub [Polynomial<8192>; 4],
    ct0is: pub [Polynomial<8192>; 4],
    ct1is: pub [Polynomial<8192>; 4],
    u: Polynomial<8192>,
    e0: Polynomial<8192>,
    e1: Polynomial<8192>,
    k1: Polynomial<8192>,
    r1is: [Polynomial<16383>; 4],
    r2is: [Polynomial<8191>; 4],
    p1is: [Polynomial<16383>; 4],
    p2is: [Polynomial<8191>; 4],
) {
    let circuit = Greco::new(
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
    circuit.verify_correct_ciphertext_encryption();
}
