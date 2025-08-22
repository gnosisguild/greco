use crate::math::polynomial::Polynomial;
use super::Params;
use safe::safe::SafeSponge;

/// BfvPkEncryption circuit structure with flattened inputs
pub struct BfvPkEncryptionCircuit<let N: u32, let L: u32> {
    params: Params<N, L>,
    // Flattened inputs instead of arrays of polynomials
    pk0is: [Field; N * L],
    pk1is: [Field; N * L],
    ct0is: [Field; N * L],
    ct1is: [Field; N * L],
    u: Polynomial<N>,
    e0: Polynomial<N>,
    e1: Polynomial<N>,
    k1: Polynomial<N>,
    r1is: [Field; (2 * N - 1) * L],
    r2is: [Field; (N - 1) * L],
    p1is: [Field; (2 * N - 1) * L],
    p2is: [Field; (N - 1) * L],
}

impl<let N: u32, let L: u32> BfvPkEncryptionCircuit<N, L> {
    /// Creates a new circuit with flattened inputs
    pub fn new(
        params: Params<N, L>,
        pk0is: [Field; N * L],
        pk1is: [Field; N * L],
        ct0is: [Field; N * L],
        ct1is: [Field; N * L],
        u: Polynomial<N>,
        e0: Polynomial<N>,
        e1: Polynomial<N>,
        k1: Polynomial<N>,
        r1is: [Field; (2 * N - 1) * L],
        r2is: [Field; (N - 1) * L],
        p1is: [Field; (2 * N - 1) * L],
        p2is: [Field; (N - 1) * L],
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

    /// Helper function to extract a polynomial from flattened array
    fn extract_polynomial(flattened: [Field; N * L], index: u32) -> [Field; N] {
        let mut coeffs = [0; N];
        let start = index * N;
        for j in 0..N {
            coeffs[j] = flattened[start + j];
        }
        coeffs
    }

    /// Get pk0is polynomial at index i
    fn get_pk0i(self, i: u32) -> Polynomial<N> {
        let coeffs = Self::extract_polynomial(self.pk0is, i);
        Polynomial::new(coeffs)
    }

    /// Get pk1is polynomial at index i
    fn get_pk1i(self, i: u32) -> Polynomial<N> {
        let coeffs = Self::extract_polynomial(self.pk1is, i);
        Polynomial::new(coeffs)
    }

    /// Get ct0is polynomial at index i
    fn get_ct0i(self, i: u32) -> Polynomial<N> {
        let coeffs = Self::extract_polynomial(self.ct0is, i);
        Polynomial::new(coeffs)
    }

    /// Get ct1is polynomial at index i
    fn get_ct1i(self, i: u32) -> Polynomial<N> {
        let coeffs = Self::extract_polynomial(self.ct1is, i);
        Polynomial::new(coeffs)
    }

    /// Get r1is polynomial at index i
    fn get_r1i(self, i: u32) -> Polynomial<(2 * N) - 1> {
        let mut coeffs = [0; (2 * N) - 1];
        let start = i * ((2 * N) - 1);
        for j in 0..((2 * N) - 1) {
            coeffs[j] = self.r1is[start + j];
        }
        Polynomial::new(coeffs)
    }

    /// Get r2is polynomial at index i
    fn get_r2i(self, i: u32) -> Polynomial<N - 1> {
        let mut coeffs = [0; N - 1];
        let start = i * (N - 1);
        for j in 0..(N - 1) {
            coeffs[j] = self.r2is[start + j];
        }
        Polynomial::new(coeffs)
    }

    /// Get p1is polynomial at index i
    fn get_p1i(self, i: u32) -> Polynomial<(2 * N) - 1> {
        let mut coeffs = [0; (2 * N) - 1];
        let start = i * ((2 * N) - 1);
        for j in 0..((2 * N) - 1) {
            coeffs[j] = self.p1is[start + j];
        }
        Polynomial::new(coeffs)
    }

    /// Get p2is polynomial at index i
    fn get_p2i(self, i: u32) -> Polynomial<N - 1> {
        let mut coeffs = [0; N - 1];
        let start = i * (N - 1);
        for j in 0..(N - 1) {
            coeffs[j] = self.p2is[start + j];
        }
        Polynomial::new(coeffs)
    }

    /// Assigns all the witness polynomial coefficients to an array (now using flattened inputs directly)
    fn payload<let SIZE: u32>(self) -> [Field; SIZE] {
        let mut inputs = [0; SIZE];
        let mut offset = 0;

        // Copy flattened arrays directly - no flattening needed
        for i in 0..(N * L) {
            inputs[offset] = self.pk0is[i];
            offset += 1;
        }
        for i in 0..(N * L) {
            inputs[offset] = self.pk1is[i];
            offset += 1;
        }
        for i in 0..(N * L) {
            inputs[offset] = self.ct0is[i];
            offset += 1;
        }
        for i in 0..(N * L) {
            inputs[offset] = self.ct1is[i];
            offset += 1;
        }
        for i in 0..N {
            inputs[offset] = self.e0.coefficients[i];
            offset += 1;
        }
        for i in 0..N {
            inputs[offset] = self.e1.coefficients[i];
            offset += 1;
        }
        for i in 0..N {
            inputs[offset] = self.u.coefficients[i];
            offset += 1;
        }
        for i in 0..N {
            inputs[offset] = self.k1.coefficients[i];
            offset += 1;
        }
        for i in 0..((2 * N - 1) * L) {
            inputs[offset] = self.r1is[i];
            offset += 1;
        }
        for i in 0..((N - 1) * L) {
            inputs[offset] = self.r2is[i];
            offset += 1;
        }
        for i in 0..((2 * N - 1) * L) {
            inputs[offset] = self.p1is[i];
            offset += 1;
        }
        for i in 0..((N - 1) * L) {
            inputs[offset] = self.p2is[i];
            offset += 1;
        }

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
            let pk0i = self.get_pk0i(i);
            let pk1i = self.get_pk1i(i);
            let r1i = self.get_r1i(i);
            let r2i = self.get_r2i(i);
            let p1i = self.get_p1i(i);
            let p2i = self.get_p2i(i);

            pk0i.range_check_1bound(self.params.pk_bounds[i]);
            pk1i.range_check_1bound(self.params.pk_bounds[i]);
            r1i.range_check_2bounds(self.params.r1_up_bounds[i], self.params.r1_low_bounds[i]);
            r2i.range_check_1bound(self.params.r2_bounds[i]);
            p1i.range_check_1bound(self.params.p1_bounds[i]);
            p2i.range_check_1bound(self.params.p2_bounds[i]);
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

            // Get polynomials and evaluate at gamma
            let pk0i = self.get_pk0i(i);
            let r1i = self.get_r1i(i);
            let r2i = self.get_r2i(i);

            let pk0is_at_gamma = pk0i.eval(gamma);
            let r1i_at_gamma = r1i.eval(gamma);
            let r2i_at_gamma = r2i.eval(gamma);

            // First step
            let pk0_u = (pk0is_at_gamma * u_at_gamma) + e0_at_gamma;

            // rhs = pk0_u + k1(gamma) * k0i
            let rhs = pk0_u + (k1_at_gamma * self.params.k0is[i]);

            // rhs = rhs + r1i(gamma) * qi
            let rhs = rhs + (r1i_at_gamma * self.params.qis[i]);

            // rhs = rhs + r2i(gamma) * cyclo(gamma)
            let rhs = rhs + r2i_at_gamma * cyclo_at_gamma;

            let ct0i = self.get_ct0i(i);
            let lhs = ct0i.eval(gamma);

            // LHS(gamma) = RHS(gamma)
            assert_eq(lhs, rhs);

            // Gamma evaluation
            let gamma = gammas.get(i + L);

            // cyclo poly is equal to x^N + 1
            let cyclo_at_gamma = gamma.pow_32(N as Field) + 1;

            let u_at_gamma = self.u.eval(gamma);
            let e1_at_gamma = self.e1.eval(gamma);

            // Get polynomials and evaluate at gamma
            let pk1i = self.get_pk1i(i);
            let p1i = self.get_p1i(i);
            let p2i = self.get_p2i(i);

            let pk1is_at_gamma = pk1i.eval(gamma);
            let p1is_at_gamma = p1i.eval(gamma);
            let p2is_at_gamma = p2i.eval(gamma);

            // Second step
            let pk1_u = (pk1is_at_gamma * u_at_gamma) + e1_at_gamma;

            //rhs = pk1_u + p2i * cyclo(gamma)
            let rhs = pk1_u + p2is_at_gamma * cyclo_at_gamma;

            //rhs = rhs + p1s * qi
            let rhs = rhs + (p1is_at_gamma * self.params.qis[i]);

            let ct1i = self.get_ct1i(i);
            let lhs = ct1i.eval(gamma);

            // LHS(gamma) = RHS(gamma)
            assert_eq(lhs, rhs);
        }
        safe.finish();
    }
}
