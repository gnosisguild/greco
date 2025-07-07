pub mod pk_encryption;

/// BfvPkEncryption params
pub struct Params<let N: u32, let L: u32> {
    q_mod_t: Field,
    pk_bounds: [u64; L],
    e_bound: u64,
    u_bound: u64,
    r1_low_bounds: [i64; L],
    r1_up_bounds: [u64; L],
    r2_bounds: [u64; L],
    p1_bounds: [u64; L],
    p2_bounds: [u64; L],
    k1_low_bound: i64,
    k1_up_bound: u64,
    qis: [Field; L],
    k0is: [Field; L],
    size: u32,
    tag: Field,
}

impl<let N: u32, let L: u32> Params<N, L> {
    pub fn new(
        q_mod_t: Field,
        pk_bounds: [u64; L],
        e_bound: u64,
        u_bound: u64,
        r1_low_bounds: [i64; L],
        r1_up_bounds: [u64; L],
        r2_bounds: [u64; L],
        p1_bounds: [u64; L],
        p2_bounds: [u64; L],
        k1_low_bound: i64,
        k1_up_bound: u64,
        qis: [Field; L],
        k0is: [Field; L],
        size: u32,
        tag: Field,
    ) -> Self {
        Params {
            q_mod_t,
            pk_bounds,
            e_bound,
            u_bound,
            r1_low_bounds,
            r1_up_bounds,
            r2_bounds,
            p1_bounds,
            p2_bounds,
            k1_low_bound,
            k1_up_bound,
            qis,
            k0is,
            size,
            tag,
        }
    }
}
