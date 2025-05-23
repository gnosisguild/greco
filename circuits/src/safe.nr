use super::pk_enc_constants_1024_2x52_2048::TAG;
use poseidon::poseidon2_permutation;

global RATE: u32 = 3;
global CAPACITY: u32 = 1;
global WIDTH: u32 = 4;

pub struct SafeSponge<let L: u32, let S: u32> {
    state: [Field; WIDTH],
    out: Vec<Field>,
    absorb_pos: u32,
    squeeze_pos: u32,
    io_pattern: [u32; L],
    io_count: u32,
}

impl<let L: u32, let S: u32> SafeSponge<L, S> {
    pub fn start(pattern: [u32; L]) -> SafeSponge<L, S> {
        let mut sponge = SafeSponge::<L, S> {
            state: [0; WIDTH],
            out: Vec::new(),
            absorb_pos: 0,
            squeeze_pos: 0,
            io_pattern: pattern,
            io_count: 0,
        };
        sponge.state[0] = TAG;
        sponge
    }

    pub fn absorb(mut self, input: [Field; S]) -> SafeSponge<L, S> {
        assert(self.io_pattern[self.io_count] as u32 == S);

        for i in 0..self.io_pattern[self.io_count] {
            if self.absorb_pos == RATE {
                self.state = poseidon2_permutation(self.state, self.state.len());
                self.absorb_pos = 0;
            }
            let pos = self.absorb_pos + CAPACITY;
            self.state[pos] = self.state[pos] + input[i];
            self.absorb_pos += 1;
        }
        self.io_count += 1;
        self.squeeze_pos = RATE;
        self
    }

    pub fn squeeze(mut self) -> Vec<Field> {
        for _ in 0..self.io_pattern[self.io_count] {
            if self.squeeze_pos == RATE {
                self.state = poseidon2_permutation(self.state, self.state.len());
                self.squeeze_pos = 0;
                self.absorb_pos = 0;
            }
            self.out.push(self.state[self.squeeze_pos + CAPACITY]);
            self.squeeze_pos += 1;
        }
        self.io_count += 1;
        self.out
    }

    pub fn finish(mut self) {
        // Clear the state
        self.state = [0; WIDTH];
        self.out = Vec::new();
        self.io_count = 0;
        self.io_pattern = [0; L];
        self.squeeze_pos = 0;
        self.absorb_pos = 0;
    }
}

#[test]
fn test_safe() {
    let pattern = [5, 4];
    let mut safe = SafeSponge::start(pattern);
    safe = safe.absorb([1, 2, 3, 4, 5]);
    let res = safe.squeeze();
    safe.finish();
    assert_eq(res.get(0), 0x0ef061c5b88d5d81e9c7c306701e01e03c5040a3d3c4b051ada3f150b44dc595);
}
