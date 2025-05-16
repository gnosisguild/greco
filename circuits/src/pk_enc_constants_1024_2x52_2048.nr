/// `N` is the degree of the cyclotomic polynomial defining the ring `Rq = Zq[X]/(X^N + 1)`.
pub global N: u32 = 1024;
/// `L` is the dimension size of the polynomials. 
pub global L: u32 = 2;
/// The coefficients of the polynomial `pk0is` and `pk1is` should exist in the interval `[-PK_BOUND, PK_BOUND]`.
pub global PK_BOUND: [u64; 2] = [2251799812767744, 2251799813160960];
/// The coefficients of the polynomial `e` should exist in the interval `[-E_BOUND, E_BOUND]` where `E_BOUND` is the upper bound of the gaussian distribution with Sigma = 3.2.
pub global E_BOUND: u64 = 19;
/// The coefficients of the polynomial `u` should exist in the interval `[-S_BOUND, S_BOUND]`.
pub global U_BOUND: u64 = 19;
/// The coefficients of the polynomials `r1is` should exist in the interval `[R1_LOW_BOUNDS[i], R1_UP_BOUNDS[i]]` where R1_LOW_BOUNDS is equal to $\frac{\frac{-(t - 1)}{2} \cdot |K_{0,i}| - (N \cdot B +2) \cdot \frac{q_i - 1}{2} + B}{q_i}` and `R1_UP_BOUNDS[i]` is equal to `$\frac{\frac{(t - 1)}{2} \cdot |K_{0,i}| + (N \cdot +2) \cdot \frac{q_i - 1}{2} + B}{q_i}` .
pub global R1_LOW_BOUNDS: [i64; 2] = [-9729, -9729];
pub global R1_UP_BOUNDS: [u64; 2] = [9729, 9729];
/// The coefficients of the polynomials `r2is` should exist in the interval `[-R2_BOUND[i], R2_BOUND[i]]` where `R2_BOUND[i]` is equal to `(qi-1)/2`.
pub global R2_BOUNDS: [u64; 2] = [2251799812767744, 2251799813160960];
/// The coefficients of the polynomials `p1is` should exist in the interval `[-P1_BOUND[i], P1_BOUND[i]]` where `P1_BOUND[i]` is equal to (((qis[i] - 1) / 2) * (N \cdot B + 2) + B ) / qis[i].
pub global P1_BOUNDS: [u64; 2] = [9728, 9728];
/// The coefficients of the polynomials `p2is` should exist in the interval `[-P2_BOUND[i], P2_BOUND[i]]` where `P2_BOUND[i]` is equal to (qis[i] - 1) / 2.
pub global P2_BOUNDS: [u64; 2] = [2251799812767744, 2251799813160960];
/// The coefficients of `k1` should exist in the interval `[K1_LOW_BOUND, K1_UP_BOUND]` where `K1_LOW_BOUND` is equal to `(-(t-1))/2` and K1_UP_BOUND` is equal to `(t-1)/2`.
pub global K1_LOW_BOUND: i64 = -1024;
pub global K1_UP_BOUND: u64 = 1023;
/// List of scalars `qis` such that `qis[i]` is the modulus of the i-th CRT basis of `q` (ciphertext space modulus).
pub global QIS: [Field; 2] = [4503599625535489, 4503599626321921];
/// List of scalars `k0is` such that `k0i[i]` is equal to the negative of the multiplicative inverses of t mod qi.
pub global K0IS: [Field; 2] = [2199023254656, 2199023255040];
