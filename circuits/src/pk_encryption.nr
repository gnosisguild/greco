use super::pk_enc_constants_1024_2x52_2048::{
    E_BOUND, K0IS, K1_LOW_BOUND, K1_UP_BOUND, L, N, P1_BOUNDS, P2_BOUNDS, PK_BOUND, QIS,
    R1_LOW_BOUNDS, R1_UP_BOUNDS, R2_BOUNDS, U_BOUND,
};
use super::polynomial::Polynomial;
use poseidon::poseidon2::Poseidon2;

pub struct BfvPkEncryptionCircuit {
    pub pk0i: [Polynomial<N>; L], // 1025, 1025
    pub pk1i: [Polynomial<N>; L], // 1025, 1025
    u: Polynomial<N>, // 1025
    e0: Polynomial<N>, // 1025
    e1: Polynomial<N>, // 1025
    k1: Polynomial<N>, // 1025
    r2is: [Polynomial<N - 1>; L], // 1024, 1024
    r1is: [Polynomial<2 * N - 1>; L], // 2048, 2048
    p2is: [Polynomial<N - 1>; L], // 1024, 1024
    p1is: [Polynomial<2 * N - 1>; L], // 2048, 2048
    pub ct0is: [Polynomial<N>; L], // 1025, 1025
    pub ct1is: [Polynomial<N>; L], // 1025, 1025
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

        for i in 0..L {
            self.pk0i[i].range_check_1bound(PK_BOUND[i]);
            self.pk1i[i].range_check_1bound(PK_BOUND[i]);
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
        let pk0i_at_gamma = self.pk0i.map(|pk| pk.eval(gamma));
        let pk1i_at_gamma = self.pk1i.map(|pk| pk.eval(gamma));
        let r1i_at_gamma = self.r1is.map(|r1| r1.eval(gamma));
        let r2i_at_gamma = self.r2is.map(|r2| r2.eval(gamma));
        let p1is_at_gamma = self.p1is.map(|p1| p1.eval(gamma));
        let p2is_at_gamma = self.p2is.map(|p2| p2.eval(gamma));

        // For each `i` Prove that LHS(gamma) = RHS(gamma)
        // pk0_u = pk0i(gamma) * u(gamma) + e0(gamma)
        // LHS = ct0i(gamma)
        // RHS = pk0_u  + k1(gamma) * k0i + r1i(gamma) * qi + r2i(gamma) * cyclo(gamma)

        // CORRECT ENCRYPTION CONSTRAINT
        for i in 0..L {
            // First step
            let pk0_u = (pk0i_at_gamma[i] * u_at_gamma) + e0_at_gamma;

            // rhs = pk0_u + k1(gamma) * k0i
            let rhs = pk0_u + (k1_at_gamma * K0IS[i]);

            // rhs = rhs + r1i(gamma) * qi
            let rhs = rhs + (r1i_at_gamma[i] * QIS[i]);

            // rhs = rhs + r2i(gamma) * cyclo(gamma) TODO: fix dummy cyclo
            let rhs = rhs + r2i_at_gamma[i];
            let lhs = self.ct0is[i].eval(gamma);

            // LHS(gamma) = RHS(gamma)
            assert_eq(lhs, rhs);

            // Second step
            let pk1_u = (pk1i_at_gamma[i] * u_at_gamma) + e1_at_gamma;

            //rhs = pk1_u + p2i * cyclo(gamma) TODO: Fix dummy cyclo
            let rhs = pk1_u + p2is_at_gamma[i];

            //rhs = rhs + p1s * qi
            let rhs = rhs + (p1is_at_gamma[i] * QIS[i]);
            let lhs = self.ct1is[i].eval(gamma);

            // LHS(gamma) = RHS(gamma)
            assert_eq(lhs, rhs);
        }
    }
}
