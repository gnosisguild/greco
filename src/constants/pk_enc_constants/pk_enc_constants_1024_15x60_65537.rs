/// `N` is the degree of the cyclotomic polynomial defining the ring `Rq = Zq[X]/(X^N + 1)`.
pub const N: usize = 1024;
/// The coefficients of the polynomial `pk0is` and `pk1is` should exist in the interval `[-PK_BOUND, PK_BOUND]`.
pub const PK_BOUND: [u64; 15] = [
    576460752303292416,
    576460752299360256,
    576460752298508288,
    576460752297984000,
    576460752297820160,
    576460752296706048,
    576460752296411136,
    576460752296214528,
    576460752294969344,
    576460752293265408,
    576460752292773888,
    576460752291823616,
    576460752290938880,
    576460752290709504,
    576460752290447360,
];
/// The coefficients of the polynomial `e` should exist in the interval `[-E_BOUND, E_BOUND]` where `E_BOUND` is the upper bound of the gaussian distribution with 𝜎 = 3.2.
pub const E_BOUND: u64 = 19;
/// The coefficients of the polynomial `u` should exist in the interval `[-S_BOUND, S_BOUND]`.
pub const U_BOUND: u64 = 19;
/// The coefficients of the polynomials `r1is` should exist in the interval `[R1_LOW_BOUNDS[i], R1_UP_BOUNDS[i]]` where R1_LOW_BOUNDS is equal to $\frac{\frac{-(t - 1)}{2} \cdot |K_{0,i}| - (N \cdot B +2) \cdot \frac{q_i - 1}{2} + B}{q_i}` and `R1_UP_BOUNDS[i]` is equal to `$\frac{\frac{(t - 1)}{2} \cdot |K_{0,i}| + (N \cdot +2) \cdot \frac{q_i - 1}{2} + B}{q_i}` .
pub const R1_LOW_BOUNDS: [i64; 15] = [
    -41688, -30870, -41317, -41479, -24000, -25422, -39592, -27470, -18559, -23996, -18968, -37075,
    -11572, -26347, -27100,
];
pub const R1_UP_BOUNDS: [u64; 15] = [
    41688, 30870, 41317, 41479, 24000, 25422, 39592, 27470, 18559, 23996, 18968, 37075, 11572,
    26347, 27100,
];
/// The coefficients of the polynomials `r2is` should exist in the interval `[-R2_BOUND[i], R2_BOUND[i]]` where `R2_BOUND[i]` is equal to `(qi-1)/2`.
pub const R2_BOUNDS: [u64; 15] = [
    576460752303292416,
    576460752299360256,
    576460752298508288,
    576460752297984000,
    576460752297820160,
    576460752296706048,
    576460752296411136,
    576460752296214528,
    576460752294969344,
    576460752293265408,
    576460752292773888,
    576460752291823616,
    576460752290938880,
    576460752290709504,
    576460752290447360,
];
/// The coefficients of the polynomials `p1is` should exist in the interval `[-P1_BOUND[i], P1_BOUND[i]]` where `P1_BOUND[i]` is equal to (((qis[i] - 1) / 2) * (N \cdot B + 2) + B ) / qis[i].
pub const P1_BOUNDS: [u64; 15] = [
    9728, 9728, 9728, 9728, 9728, 9728, 9728, 9728, 9728, 9728, 9728, 9728, 9728, 9728, 9728,
];
/// The coefficients of the polynomials `p2is` should exist in the interval `[-P2_BOUND[i], P2_BOUND[i]]` where `P2_BOUND[i]` is equal to (qis[i] - 1) / 2.
pub const P2_BOUNDS: [u64; 15] = [
    576460752303292416,
    576460752299360256,
    576460752298508288,
    576460752297984000,
    576460752297820160,
    576460752296706048,
    576460752296411136,
    576460752296214528,
    576460752294969344,
    576460752293265408,
    576460752292773888,
    576460752291823616,
    576460752290938880,
    576460752290709504,
    576460752290447360,
];
/// The coefficients of `k1` should exist in the interval `[K1_LOW_BOUND, K1_UP_BOUND]` where `K1_LOW_BOUND` is equal to `(-(t-1))/2` and K1_UP_BOUND` is equal to `(t-1)/2`.
pub const K1_LOW_BOUND: i64 = -32768;
pub const K1_UP_BOUND: u64 = 32768;
/// List of scalars `qis` such that `qis[i]` is the modulus of the i-th CRT basis of `q` (ciphertext space modulus).
pub const QIS: [&str; 15] = [
    "1152921504606584833",
    "1152921504598720513",
    "1152921504597016577",
    "1152921504595968001",
    "1152921504595640321",
    "1152921504593412097",
    "1152921504592822273",
    "1152921504592429057",
    "1152921504589938689",
    "1152921504586530817",
    "1152921504585547777",
    "1152921504583647233",
    "1152921504581877761",
    "1152921504581419009",
    "1152921504580894721",
];
/// List of scalars `k0is` such that `k0i[i]` is equal to the negative of the multiplicative inverses of t mod qi.
pub const K0IS: [&str; 15] = [
    "1124457781908666798",
    "743839052427601194",
    "1111422170948171465",
    "1117121952253736973",
    "502126104424574846",
    "552157518114552474",
    "1050730055179823439",
    "624214012656257690",
    "310690856959624444",
    "501985369079705462",
    "325081045565655086",
    "962154749991507364",
    "64878992155545191",
    "584702565692244436",
    "611213585534257079",
];
