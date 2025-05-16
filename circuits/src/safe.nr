use poseidon::poseidon::bn254::hash_2;
use std::hash::hash_to_field;
use std::hash::poseidon2_permutation;

struct SafeSponge<let L: u32> {
    state: [Field; 5],
    absorb_pos: u32,
    squeeze_pos: u32,
    // expected IO pattern (compact 32-bit words)
    io_pattern: [u32; L],
    io_count: u32,
}

impl<let L: u32> SafeSponge<L> {
    pub fn start(pattern: [u32; L], domain_sep: Field) -> SafeSponge<L> {
        let mut tag = hash_to_field(pattern.map(|x| x as Field));
        tag = tag + domain_sep;

        let mut sponge = SafeSponge::<L> {
            state: [0; 5],
            absorb_pos: 0,
            squeeze_pos: 0,
            io_pattern: pattern,
            io_count: 0,
        };
        // TODO: Maybe to index 4?
        sponge.state[0] = tag;
        sponge
    }

    pub fn absorb<let M: u32>(mut self, input: [Field; M]) {
        for i in 0..M {
            if self.absorb_pos == 4 {
                self.state = poseidon_permutation(self.state);
                self.absorb_pos = 0;
            }
            self.state[self.absorb_pos] = self.state[self.absorb_pos] + input[i];
            self.absorb_pos += 1;
        }
        // TODO: Add io check and correct this part
        assert(self.io_pattern[self.io_count] == M);
        self.io_count = self.io_count + 1;
        self.squeeze_pos = 4;
    }

    pub fn squeeze<let M: u32>(mut self) -> [Field; M] {
        let mut out = [0; M];
        for i in 0..M {
            if self.squeeze_pos == 4 {
                self.state = poseidon2_permutation(self.state, self.state.len());
                self.squeeze_pos = 0;
                self.absorb_pos = 0;
            }
            out[i] = self.state[self.squeeze_pos + 1];
            self.squeeze_pos = self.squeeze_pos + 1;
        }
        //assert(self.io_pattern[self.io_count] == M);
        self.io_count += 1;
        out
    }

    pub fn finish(mut self) {
        assert(self.io_count == self.io_pattern.len());
        // Clear the state
        self.state = [0; 5];
        self.io_count = 0;
        self.io_pattern = [0; L];
        self.squeeze_pos = 0;
        self.absorb_pos = 0;
    }
}

fn poseidon_permutation(state: [Field; 5]) -> [Field; 5] {
    let mut new_state = state;
    for i in 0..4 {
        new_state[4] = hash_2([new_state[i], new_state[4]]);
    }
    new_state
}

// #[test]
// fn test_safe() {
//     let mut sponge = SafeSponge::start([1, 2], 0);
//     sponge.absorb([1]);
//     let res: [Field; 1] = sponge.squeeze();
// }
