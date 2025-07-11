use crate::math::polynomial::Polynomial;
use super::Params;
use safe::safe::SafeSponge;

/// BfvPkEncryption circuit structure
pub struct BfvPkEncryptionCircuit<let N: u32, let L: u32> {
    params: Params<N, L>,
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

impl<let N: u32, let L: u32> BfvPkEncryptionCircuit<N, L> {
    /// Creates a new circuit
    pub fn new(
        params: Params<N, L>,
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
    ) -> BfvPkEncryptionCircuit<N, L> {
        BfvPkEncryptionCircuit {
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
        }
    }

    /// Assigns all the witness polynomial coefficients to an array
    fn payload<let SIZE: u32>(self) -> [Field; SIZE] {
        let mut inputs = [0; SIZE];
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
        let (inputs, _) = flatten(inputs, self.p2is, offset);

        inputs
    }

    pub fn correct_encryption(self) {
        // Binary check

        assert(
            0 == self.k1.coefficients[N - 1] * (self.params.q_mod_t - self.k1.coefficients[N - 1]),
        );

        // Bound check
        self.u.range_check_1bound(self.params.u_bound);
        self.e0.range_check_1bound(self.params.e_bound);
        self.e1.range_check_1bound(self.params.e_bound);
        self.k1.range_check_2bounds(self.params.k1_up_bound, self.params.k1_low_bound);

        for i in 0..L {
            self.pk0is[i].range_check_1bound(self.params.pk_bounds[i]);
            self.pk1is[i].range_check_1bound(self.params.pk_bounds[i]);
            self.r1is[i].range_check_2bounds(
                self.params.r1_up_bounds[i],
                self.params.r1_low_bounds[i],
            );
            self.r2is[i].range_check_1bound(self.params.r2_bounds[i]);
            self.p1is[i].range_check_1bound(self.params.p1_bounds[i]);
            self.p2is[i].range_check_1bound(self.params.p2_bounds[i]);
        }

        // We assign all the coefficients to an array to generate challenge values
        // SIZE formula is (10 * N - 4) * L + 4 * N
        let inputs = self.payload::<(10 * N - 4) * L + 4 * N>();
        let mut safe = SafeSponge::start([(10 * N - 4) * L + 4 * N, 2 * L], self.params.tag);
        safe = safe.absorb(inputs);
        let gammas = safe.squeeze();

        // CORRECT ENCRYPTION CONSTRAINT
        // For each `i` Prove that LHS(gamma) = RHS(gamma)
        // pk0_u = pk0is(gamma) * u(gamma) + e0(gamma)
        // LHS = ct0i(gamma)
        // RHS = pk0_u  + k1(gamma) * k0i + r1i(gamma) * qi + r2i(gamma) * cyclo(gamma)
        for i in 0..L {
            let gamma = gammas.get(i);

            // cyclo poly is equal to x^N + 1
            let cyclo_at_gamma = gamma.pow_32(N as Field) + 1;

            // Gamma evaluation
            let u_at_gamma = self.u.eval(gamma);
            let e0_at_gamma = self.e0.eval(gamma);
            let k1_at_gamma = self.k1.eval(gamma);
            let pk0is_at_gamma = self.pk0is.map(|pk| pk.eval(gamma));
            let r1i_at_gamma = self.r1is.map(|r1| r1.eval(gamma));
            let r2i_at_gamma = self.r2is.map(|r2| r2.eval(gamma));

            // First step
            let pk0_u = (pk0is_at_gamma[i] * u_at_gamma) + e0_at_gamma;

            // rhs = pk0_u + k1(gamma) * k0i
            let rhs = pk0_u + (k1_at_gamma * self.params.k0is[i]);

            // rhs = rhs + r1i(gamma) * qi
            let rhs = rhs + (r1i_at_gamma[i] * self.params.qis[i]);

            // rhs = rhs + r2i(gamma) * cyclo(gamma)
            let rhs = rhs + r2i_at_gamma[i] * cyclo_at_gamma;
            let lhs = self.ct0is[i].eval(gamma);

            // LHS(gamma) = RHS(gamma)
            assert_eq(lhs, rhs);

            // Gamma evaluation
            let gamma = gammas.get(i + L);

            // cyclo poly is equal to x^N + 1
            let cyclo_at_gamma = gamma.pow_32(N as Field) + 1;

            let u_at_gamma = self.u.eval(gamma);
            let e1_at_gamma = self.e1.eval(gamma);
            let pk1is_at_gamma = self.pk1is.map(|pk| pk.eval(gamma));
            let p1is_at_gamma = self.p1is.map(|p1| p1.eval(gamma));
            let p2is_at_gamma = self.p2is.map(|p2| p2.eval(gamma));

            // Second step
            let pk1_u = (pk1is_at_gamma[i] * u_at_gamma) + e1_at_gamma;

            //rhs = pk1_u + p2i * cyclo(gamma)
            let rhs = pk1_u + p2is_at_gamma[i] * cyclo_at_gamma;

            //rhs = rhs + p1s * qi
            let rhs = rhs + (p1is_at_gamma[i] * self.params.qis[i]);
            let lhs = self.ct1is[i].eval(gamma);

            // LHS(gamma) = RHS(gamma)
            assert_eq(lhs, rhs);
        }
        safe.finish();
    }
}

/// Flattens an array of polynomials into a 1D array of Field elements.
///
/// # Generic Parameters
/// * `A` - The size of each polynomial in the array
/// * `L` - The number of polynomials in the array
///
/// # Arguments
/// * `inputs` - The target array to store flattened coefficients
/// * `poly` - Array of polynomials to flatten
/// * `offset` - Starting position in the inputs array
///
/// # Returns
/// A tuple containing:
/// * The updated inputs array with flattened coefficients
/// * The new offset position
/// ```
fn flatten<let A: u32, let L: u32, let SIZE: u32>(
    mut inputs: [Field; SIZE],
    poly: [Polynomial<A>; L],
    mut offset: u32,
) -> ([Field; SIZE], u32) {
    for j in 0..L {
        for i in 0..A {
            inputs[i + offset] = poly[j].coefficients[i];
        }
        offset += A;
    }

    (inputs, offset)
}
