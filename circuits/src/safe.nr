use super::pk_enc_constants_1024_2x52_2048::TAG;
use poseidon::poseidon2_permutation;

global RATE: u32 = 2;
global CAPACITY: u32 = 1;
global WIDTH: u32 = 4;

struct SafeSponge<let L: u32, let S: u32> {
    state: [Field; WIDTH],
    absorb_pos: u32,
    squeeze_pos: u32,
    io_pattern: [i32; L],
    io_count: u32,
}

impl<let L: u32, let S: u32> SafeSponge<L, S> {
    pub fn start(pattern: [i32; L]) -> SafeSponge<L, S> {
        let mut sponge = SafeSponge::<L, S> {
            state: [0; WIDTH],
            absorb_pos: 0,
            squeeze_pos: 0,
            io_pattern: pattern,
            io_count: 0,
        };
        sponge.state[0] = TAG;
        sponge
    }

    pub fn absorb(mut self, input: [Field; S]) -> SafeSponge<L, S> {
        assert(self.io_pattern[self.io_count] > 0);
        assert(self.io_pattern[self.io_count] as u32 == input.len());

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

    pub fn squeeze(mut self) -> Field {
        assert(self.io_pattern[self.io_count] < 0);
        let mut out = 0;

        for _ in 0..-(self.io_pattern[self.io_count]) {
            if self.squeeze_pos == RATE {
                self.state = poseidon2_permutation(self.state, self.state.len());
                self.squeeze_pos = 0;
                self.absorb_pos = 0;
            }
            self.squeeze_pos += 1;
            out = out + self.state[self.squeeze_pos + CAPACITY];
        }
        self.io_count += 1;
        out
    }

    pub fn finish(mut self) {
        assert(self.io_count == self.io_pattern.len());
        // Clear the state
        self.state = [0; WIDTH];
        self.io_count = 0;
        self.io_pattern = [0; L];
        self.squeeze_pos = 0;
        self.absorb_pos = 0;
    }
}

#[test]
fn test_safe() {
    let pattern = [5, -1];
    let mut safe = SafeSponge::start(pattern);
    safe = safe.absorb([1, 2, 3, 4, 5]);
    let res = safe.squeeze();
    assert_eq(res, 0x261007ad6fed700ec20bfc988556ed9e9fbf84fdca7eaeec3eb644100ffbfb83);
}
