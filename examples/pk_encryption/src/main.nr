use circuits::constants::{L, N};
use circuits::crypto::{Params, pk_encryption::BfvPkEncryptionCircuit};
use circuits::math::polynomial::Polynomial;

fn main(
    params: pub Params<N, L>,
    pk0is: pub [Polynomial<N>; L],
    pk1is: pub [Polynomial<N>; L],
    ct0is: pub [Polynomial<N>; L],
    ct1is: pub [Polynomial<N>; L],
    u: Polynomial<N>,
    e0: Polynomial<N>,
    e1: Polynomial<N>,
    k1: Polynomial<N>,
    r1is: [Polynomial<(2 * N) - 1>; L],
    r2is: [Polynomial<N - 1>; L],
    p1is: [Polynomial<(2 * N) - 1>; L],
    p2is: [Polynomial<N - 1>; L],
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
