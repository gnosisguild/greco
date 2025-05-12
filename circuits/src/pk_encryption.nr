use super::pk_enc_constants_1024_2x52_2048::{
    E_BOUND, K0IS, K1_LOW_BOUND, K1_UP_BOUND, L, N, P1_BOUNDS, P2_BOUNDS, PK_BOUND, QIS,
    R1_LOW_BOUNDS, R1_UP_BOUNDS, R2_BOUNDS, U_BOUND,
};
use super::polynomial::Polynomial;
use poseidon::poseidon2::Poseidon2;

// Total size of the all polynomial's coefficients.
global size: u32 = (10 * L + 4) * N - 8;

/// BfvPkEncryption circuit structure
pub struct BfvPkEncryptionCircuit {
    pk0is: [Polynomial<N>; L],
    pk1is: [Polynomial<N>; L],
    ct0is: [Polynomial<N>; L],
    ct1is: [Polynomial<N>; L],
    u: Polynomial<N>,
    e0: Polynomial<N>,
    e1: Polynomial<N>,
    k1: Polynomial<N>,
    r1is: [Polynomial<(2 * N) - 1>; L],
    r2is: [Polynomial<N - 1>; L],
    p1is: [Polynomial<(2 * N) - 1>; L],
    p2is: [Polynomial<N - 1>; L],
}

impl BfvPkEncryptionCircuit {
    /// Creates a new circuit
    pub fn new(
        pk0is: [Polynomial<N>; L],
        pk1is: [Polynomial<N>; L],
        ct0is: [Polynomial<N>; L],
        ct1is: [Polynomial<N>; L],
        u: Polynomial<N>,
        e0: Polynomial<N>,
        e1: Polynomial<N>,
        k1: Polynomial<N>,
        r1is: [Polynomial<2 * N - 1>; L],
        r2is: [Polynomial<N - 1>; L],
        p1is: [Polynomial<2 * N - 1>; L],
        p2is: [Polynomial<N - 1>; L],
    ) -> BfvPkEncryptionCircuit {
        BfvPkEncryptionCircuit { pk0is, pk1is, ct0is, ct1is, u, e0, e1, k1, r1is, r2is, p1is, p2is }
    }

    /// Assigns all the witness polynomial coefficients to an array
    fn payload(self) -> [Field; size] {
        let mut inputs = [0; size];
        let mut offset = 0;

        let (inputs, offset) = flatten(inputs, self.pk0is, offset);
        let (inputs, offset) = flatten(inputs, self.pk1is, offset);
        let (inputs, offset) = flatten(inputs, self.ct0is, offset);
        let (inputs, offset) = flatten(inputs, self.ct1is, offset);
        let (inputs, offset) = flatten(inputs, [self.e0], offset);
        let (inputs, offset) = flatten(inputs, [self.e1], offset);
        let (inputs, offset) = flatten(inputs, [self.u], offset);
        let (inputs, offset) = flatten(inputs, [self.k1], offset);
        let (inputs, offset) = flatten(inputs, self.r1is, offset);
        let (inputs, offset) = flatten(inputs, self.r2is, offset);
        let (inputs, offset) = flatten(inputs, self.p1is, offset);
        let (inputs, _offset) = flatten(inputs, self.p2is, offset);

        inputs
    }

    pub fn correct_encryption(self) {
        // TODO: Add domain separation
        // We assign all the coefficients to an array to generate a challenge value with Poseidon hash function
        let inputs = self.payload();
        let gamma = Poseidon2::hash(inputs, size);

        // cyclo poly is equal to x^N + 1
        let cyclo_at_gamma = gamma.pow_32(N as Field) + 1;

        // Bound checks
        self.u.range_check_1bound(U_BOUND);
        self.e0.range_check_1bound(E_BOUND);
        self.e1.range_check_1bound(E_BOUND);
        self.k1.range_check_2bounds(K1_UP_BOUND, K1_LOW_BOUND);

        for i in 0..L {
            self.pk0is[i].range_check_1bound(PK_BOUND[i]);
            self.pk1is[i].range_check_1bound(PK_BOUND[i]);
            self.r1is[i].range_check_2bounds(R1_UP_BOUNDS[i], R1_LOW_BOUNDS[i]);
            self.r2is[i].range_check_1bound(R2_BOUNDS[i]);
            self.p1is[i].range_check_1bound(P1_BOUNDS[i]);
            self.p2is[i].range_check_1bound(P2_BOUNDS[i]);
        }

        // Gamma evaluation
        let u_at_gamma = self.u.eval(gamma);
        let e0_at_gamma = self.e0.eval(gamma);
        let e1_at_gamma = self.e1.eval(gamma);
        let k1_at_gamma = self.k1.eval(gamma);
        let pk0is_at_gamma = self.pk0is.map(|pk| pk.eval(gamma));
        let pk1is_at_gamma = self.pk1is.map(|pk| pk.eval(gamma));
        let r1i_at_gamma = self.r1is.map(|r1| r1.eval(gamma));
        let r2i_at_gamma = self.r2is.map(|r2| r2.eval(gamma));
        let p1is_at_gamma = self.p1is.map(|p1| p1.eval(gamma));
        let p2is_at_gamma = self.p2is.map(|p2| p2.eval(gamma));

        // For each `i` Prove that LHS(gamma) = RHS(gamma)
        // pk0_u = pk0is(gamma) * u(gamma) + e0(gamma)
        // LHS = ct0i(gamma)
        // RHS = pk0_u  + k1(gamma) * k0i + r1i(gamma) * qi + r2i(gamma) * cyclo(gamma)

        // CORRECT ENCRYPTION CONSTRAINT
        for i in 0..L {
            // First step
            let pk0_u = (pk0is_at_gamma[i] * u_at_gamma) + e0_at_gamma;

            // rhs = pk0_u + k1(gamma) * k0i
            let rhs = pk0_u + (k1_at_gamma * K0IS[i]);

            // rhs = rhs + r1i(gamma) * qi
            let rhs = rhs + (r1i_at_gamma[i] * QIS[i]);

            // rhs = rhs + r2i(gamma) * cyclo(gamma)
            let rhs = rhs + r2i_at_gamma[i] * cyclo_at_gamma;
            let lhs = self.ct0is[i].eval(gamma);

            // LHS(gamma) = RHS(gamma)
            assert_eq(lhs, rhs);

            // Second step
            let pk1_u = (pk1is_at_gamma[i] * u_at_gamma) + e1_at_gamma;

            //rhs = pk1_u + p2i * cyclo(gamma)
            let rhs = pk1_u + p2is_at_gamma[i] * cyclo_at_gamma;

            //rhs = rhs + p1s * qi
            let rhs = rhs + (p1is_at_gamma[i] * QIS[i]);
            let lhs = self.ct1is[i].eval(gamma);

            // LHS(gamma) = RHS(gamma)
            assert_eq(lhs, rhs);
        }
    }
}

/// Flattens a 2D array of `Polynomial<N>` into a 1D `inputs` array of `Field` elements.
/// # Example
/// ```rust
/// let (flattened, new_offset) = flatten(inputs, polynomials, offset);
/// ```
fn flatten<let N: u32, let L: u32>(
    mut inputs: [Field; size],
    poly: [Polynomial<N>; L],
    mut offset: u32,
) -> ([Field; size], u32) {
    for j in 0..L {
        for i in 0..N {
            inputs[i + offset] = poly[j].coefficients[i];
        }
        offset += N;
    }

    (inputs, offset)
}

#[test]
fn test_correct_encryption() {
    // TODO: Add test data
    // let circuit = BfvPkEncryptionCircuit::new(
    //     pk0is,
    //     pk1is,
    //     ct0is,
    //     ct1is,
    //     u,
    //     e0,
    //     e1,
    //     k1,
    //     r1is,
    //     r2is,
    //     p1is,
    //     p2is
    // );
    // circuit.correct_encryption();
}
