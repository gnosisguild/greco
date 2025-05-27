use super::pk_enc_constants_1024_2x52_2048::{
    E_BOUND, K0IS, K1_LOW_BOUND, K1_UP_BOUND, L, N, P1_BOUNDS, P2_BOUNDS, PK_BOUND, QIS,
    R1_LOW_BOUNDS, R1_UP_BOUNDS, R2_BOUNDS, U_BOUND,
};
use super::polynomial::Polynomial;
use super::safe::SafeSponge;

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
        // Bound check
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

        // We assign all the coefficients to an array to generate challenge values
        let inputs = self.payload();
        let mut safe = SafeSponge::start([size, 2 * L]);
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
            let rhs = pk0_u + (k1_at_gamma * K0IS[i]);

            // rhs = rhs + r1i(gamma) * qi
            let rhs = rhs + (r1i_at_gamma[i] * QIS[i]);

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
            let rhs = rhs + (p1is_at_gamma[i] * QIS[i]);
            let lhs = self.ct1is[i].eval(gamma);

            // LHS(gamma) = RHS(gamma)
            assert_eq(lhs, rhs);
        }
        safe.finish();
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

// #[test]
// // Test case for correct ecnryption function with small Polynomials
// fn correct_enc_test() {
//     let pk0is = [
//         Polynomial::new([
//             -1403802783733623,
//             -1805630581469767,
//             1074396299991494,
//             -1191821857299158,
//             232760576561784,
//             2009834201656190,
//             -1712673780202959,
//             222679397902536,
//         ]),
//         Polynomial::new([
//             -823720773323418,
//             1891450372563020,
//             -128745212030841,
//             -441455275386824,
//             -872158998692394,
//             -475552226335370,
//             -90747098345782,
//             -1036800016830931,
//         ]),
//     ];
//     let pk1is = [
//         Polynomial::new([
//             1648718294067155,
//             2245512895797282,
//             -447800124309729,
//             -2111249083329095,
//             2066133022760530,
//             -584869862223593,
//             -842704615493033,
//             142108104908355,
//         ]),
//         Polynomial::new([
//             -291899183447317,
//             795790016923348,
//             1860195844307355,
//             668408031295103,
//             -838693223974419,
//             -625968156329695,
//             608259791953513,
//             1437744972024037,
//         ]),
//     ];
//     let ct0is = [
//         Polynomial::new([
//             -370220263333386,
//             -327862368600833,
//             1884754970694558,
//             2066651122784285,
//             2013216223322597,
//             1229154658621711,
//             1045585492319234,
//             -1489882202557254,
//         ]),
//         Polynomial::new([
//             842538933607836,
//             784279925552216,
//             1527556023823927,
//             -2059593818569710,
//             1102382499342986,
//             -483217896470498,
//             -1750257971413370,
//             -893468153257829,
//         ]),
//     ];
//     let ct1is = [
//         Polynomial::new([
//             776823628252784,
//             2055137419154754,
//             2140462596874717,
//             160113978504870,
//             -135877097553202,
//             1638485724211152,
//             1543584081500201,
//             -2184470657222490,
//         ]),
//         Polynomial::new([
//             -69539390612689,
//             137155248672618,
//             2007797448277144,
//             -1336187492633734,
//             529715659543507,
//             1299911737964247,
//             -2182441341470478,
//             2086132173908181,
//         ]),
//     ];
//     let u = Polynomial::new([3, -4, -5, 3, 0, 1, -1, -3]);
//     let e0 = Polynomial::new([-1, 1, -4, -5, -3, -1, -2, 7]);
//     let e1 = Polynomial::new([-3, -2, 5, -2, 3, 1, 1, -4]);
//     let k1 = Polynomial::new([-1, 1, 1, 0, -1, 0, 0, -1]);

//     let r1is = [
//         Polynomial::new([1, 0, -4, 1, 1, -3, 4, -1, -4, 2, -1, 2, 1, -1, 0]),
//         Polynomial::new([1, -2, 1, 3, -1, -1, -2, 1, 1, -1, 0, 0, 0, -1, -1]),
//     ];
//     let r2is = [
//         Polynomial::new([
//             -292191274334620,
//             -198319390525191,
//             549673357620291,
//             -1447293409820005,
//             819704459949298,
//             633843175148384,
//             304055897534219,
//         ]),
//         Polynomial::new([
//             -2032437306351667,
//             37965041661110,
//             -670166366594408,
//             -772999718143301,
//             -963821656991458,
//             1444300343040276,
//             1625631329151772,
//         ]),
//     ];
//     let p1is = [
//         Polynomial::new([-1, 0, 4, 2, -5, 0, 4, -1, 2, 0, -1, 1, 0, 0, 0]),
//         Polynomial::new([0, -1, -1, 2, 3, -1, -3, -1, 3, 3, -1, -1, 0, 0, 1]),
//     ];
//     let p2is = [
//         Polynomial::new([
//             -442555256665976,
//             -141665511123226,
//             554644924312134,
//             1816757098462336,
//             -1100936582861016,
//             -842421660070544,
//             -1758146342497421,
//         ]),
//         Polynomial::new([
//             875697550341951,
//             948632841762609,
//             646676243856663,
//             1283007665658960,
//             -1407477911095363,
//             1080083006414097,
//             1895767463658375,
//         ]),
//     ];

//     let bfv = BfvPkEncryptionCircuit::new(
//         pk0is,
//         pk1is,
//         ct0is,
//         ct1is,
//         u,
//         e0,
//         e1,
//         k1,
//         r1is,
//         r2is,
//         p1is,
//         p2is,
//     );
//     bfv.correct_encryption();
// }
