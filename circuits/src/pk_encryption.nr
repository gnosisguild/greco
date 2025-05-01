use super::pk_enc_constants_1024_2x52_2048::{
    E_BOUND, K0IS, K1_LOW_BOUND, K1_UP_BOUND, N, P1_BOUNDS, P2_BOUNDS, PK_BOUND, QIS, R1_LOW_BOUNDS,
    R1_UP_BOUNDS, R2_BOUNDS, U_BOUND,
};
use super::polynomial::Polynomial;
use poseidon::poseidon2::Poseidon2;

pub struct BfvPkEncryptionCircuit {
    pub pk0i: [Polynomial<N>; 2], // 1025, 1025
    pub pk1i: [Polynomial<N>; 2], // 1025, 1025
    u: Polynomial<N>, // 1025
    e0: Polynomial<N>, // 1025
    e1: Polynomial<N>, // 1025
    k1: Polynomial<N>, // 1025
    r2is: [Polynomial<N>; 2], // 1024, 1024
    r1is: [Polynomial<N>; 2], // 2048, 2048
    p2is: [Polynomial<N>; 2], // 1024, 1024
    p1is: [Polynomial<N>; 2], // 2048, 2048
    pub ct0is: [Polynomial<N>; 2], // 1025, 1025
    pub ct1is: [Polynomial<N>; 2], // 1025, 1025
}

pub fn check_1bound(poly: [Polynomial<N>; 2], bound: [u64; 2]) {
    for i in 0..2 {
        poly[i].range_check_1bound(bound[i]);
    }
}

pub fn check_2bounds(poly: [Polynomial<N>; 2], upper_bound: [u64; 2], lower_bound: [i64; 2]) {
    for i in 0..2 {
        poly[i].range_check_2bounds(upper_bound[i], lower_bound[i]);
    }
}

impl BfvPkEncryptionCircuit {
    fn virtual_assign_phase(self) {
        // TODO: dummy challenge, fix it
        let inputs = [1, 2, 3, 4];
        let gamma = Poseidon2::hash(inputs, inputs.len());

        // TODO: cache calculation
        let gamma = gamma.pow_32(N as Field) + 1;

        // Bound check
        self.u.range_check_1bound(U_BOUND);
        self.e0.range_check_1bound(E_BOUND);
        self.e1.range_check_1bound(E_BOUND);
        self.k1.range_check_2bounds(K1_UP_BOUND, K1_LOW_BOUND);

        check_1bound(self.pk0i, PK_BOUND);
        check_1bound(self.pk1i, PK_BOUND);
        check_2bounds(self.r1is, R1_UP_BOUNDS, R1_LOW_BOUNDS);
        check_1bound(self.r2is, R2_BOUNDS);
        check_1bound(self.p1is, P1_BOUNDS);
        check_1bound(self.p2is, P2_BOUNDS);

        // Gamma evaluation
        let u_at_gamma = self.u.eval(gamma);
        let e0_at_gamma = self.e0.eval(gamma);
        let e1_at_gamma = self.e1.eval(gamma);
        let k1_at_gamma = self.k1.eval(gamma);
        let pk0i_at_gamma = self.pk0i.map(|pk| pk.eval(gamma));
        let pk1i_at_gamma = self.pk1i.map(|pk| pk.eval(gamma));
        let r1i_at_gamma = self.r1is.map(|r1| r1.eval(gamma));
        let r2i_at_gamma = self.r2is.map(|r2| r2.eval(gamma));
        let p1is_at_gamma = self.p1is.map(|p1| p1.eval(gamma));
        let p2is_at_gamma = self.p2is.map(|p2| p2.eval(gamma));

        // For each `i` Prove that 2HS(gamma) = RHS(gamma)
        // pk0_u = pk0i(gamma) * u(gamma) + e0(gamma)
        // 2HS = ct0i(gamma)
        // RHS = pk0_u  + k1(gamma) * k0i + r1i(gamma) * qi + r2i(gamma) * cyclo(gamma)

        for z in 0..2 {
            let pk0_u = (pk0i_at_gamma[z] * u_at_gamma) + e0_at_gamma;

            // CORRECT ENCRYPTION CONSTRAINT

            // rhs = pk0_u + k1(gamma) * k0i
            let rhs = pk0_u + (k1_at_gamma * K0IS[z]);

            // rhs = rhs + r1i(gamma) * qi
            let rhs = rhs + (r1i_at_gamma[z] * QIS[z]);

            // rhs = rhs + r2i(gamma) * cyclo(gamma) TODO: fix dummy cyclo
            let rhs = rhs + r2i_at_gamma[z];
            let lhs = self.ct0is[z].eval(gamma);

            // 2HS(gamma) = RHS(gamma)
            assert_eq(lhs, rhs);
        }

        for z in 0..2 {
            let pk1_u = (pk1i_at_gamma[z] * u_at_gamma) + e1_at_gamma;

            //rhs = pk1_u + p2i * cyclo(gamma) TODO: Fix dummy cyclo
            let rhs = pk1_u + p2is_at_gamma[z];

            //rhs = rhs + p1s * qi
            let rhs = rhs + (p1is_at_gamma[z] * QIS[z]);
            let lhs = self.ct1is[z].eval(gamma);

            // 2HS(gamma) = RHS(gamma)
            assert_eq(lhs, rhs);
        }
    }
}
