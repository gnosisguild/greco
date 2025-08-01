use crate::math::polynomial::Polynomial;
use super::Params;
use safe::safe::SafeSponge;

use std::hash::pedersen_commitment;
use std::embedded_curve_ops::EmbeddedCurvePoint;

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

	    // pk0is
	    for i in 0..L {
		let c = pedersen_commitment(self.pk0is[i].coefficients);
		let (new_inputs, new_offset) = flatten_commitment(inputs, c, offset);
		inputs = new_inputs;
		offset = new_offset;
	    }

	    println("inputs");
	    println(inputs);

	    // pk1is
	    for i in 0..L {
		let c = pedersen_commitment(self.pk1is[i].coefficients);
		let (new_inputs, new_offset) = flatten_commitment(inputs, c, offset);
		inputs = new_inputs;
		offset = new_offset;
	    }

	    // ct0is
	    for i in 0..L {
		let c = pedersen_commitment(self.ct0is[i].coefficients);
		let (new_inputs, new_offset) = flatten_commitment(inputs, c, offset);
		inputs = new_inputs;
		offset = new_offset;
	    }

	    // ct1is
	    for i in 0..L {
		let c = pedersen_commitment(self.ct1is[i].coefficients);
		let (new_inputs, new_offset) = flatten_commitment(inputs, c, offset);
		inputs = new_inputs;
		offset = new_offset;
	    }

	    // e0, e1, u, k1
	    for poly in [self.e0, self.e1, self.u, self.k1] {
		let c = pedersen_commitment(poly.coefficients);
		let (new_inputs, new_offset) = flatten_commitment(inputs, c, offset);
		inputs = new_inputs;
		offset = new_offset;
	    }

	    // r1is
	    for i in 0..L {
		let c = pedersen_commitment(self.r1is[i].coefficients);
		let (new_inputs, new_offset) = flatten_commitment(inputs, c, offset);
		inputs = new_inputs;
		offset = new_offset;
	    }

	    // r2is
	    for i in 0..L {
		let c = pedersen_commitment(self.r2is[i].coefficients);
		let (new_inputs, new_offset) = flatten_commitment(inputs, c, offset);
		inputs = new_inputs;
		offset = new_offset;
	    }

	    // p1is
	    for i in 0..L {
		let c = pedersen_commitment(self.p1is[i].coefficients);
		let (new_inputs, new_offset) = flatten_commitment(inputs, c, offset);
		inputs = new_inputs;
		offset = new_offset;
	    }

	    // p2is
	    for i in 0..L {
		let c = pedersen_commitment(self.p2is[i].coefficients);
		let (new_inputs, new_offset) = flatten_commitment(inputs, c, offset);
		inputs = new_inputs;
		offset = new_offset;
	    }
 	println("inputs");
	    println(inputs);


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
        // SIZE formula is 2 * (8 * L + 4)
        let inputs = self.payload::< 16 * L + 8 >();
        let mut safe = SafeSponge::start([ 16 * L + 8, 2 * L], self.params.tag);
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

/// Flattens Pedersen committment (which is an elliptic curve point [x,y]) into a 1D array of Field elements.
///
/// # Arguments
/// * `inputs` - The target array to store flattened coefficients
/// * `commitment` - An ellitpic curve point
/// * `offset` - Starting position in the inputs array
///
/// # Returns
/// A tuple containing:
/// * The updated inputs array with flattened values 
/// * The new offset position
/// ```
fn flatten_commitment<let SIZE: u32>(
    mut inputs: [Field; SIZE],
    commitment: EmbeddedCurvePoint,
    offset: u32,
) -> ([Field; SIZE], u32) {
    inputs[offset] = commitment.x;
    inputs[offset + 1] = commitment.y;
    (inputs, offset + 2)
}
