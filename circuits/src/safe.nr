use poseidon::poseidon::bn254::hash_2;
use std::hash::hash_to_field;
use std::hash::poseidon2_permutation;

global RATE: u32 = 4;
global CAPACITY: u32 = 1;
global WIDTH: u32 = 5;
global ABSORB_MASK: i64 = 0x80000_0000;

fn tag_input<let L: u32>(pattern: [i32; L], domain_sep: Field) -> [Field] {
    // TODO: discuss this approach, since we dont have enums negative values are squeeze and positives are absorb. This approach is not safe
    // Also we have to check that each call should be smaller than 2^31 bits, that is why we are using i32 here.
    assert(pattern[0] > 0);
    assert(pattern[L] < 0);

    let mut inputs = [0; L + 2];
    inputs[0] = ABSORB_MASK;

    for i in 0..pattern.len() {
        let mut prev = inputs[i];
        let value = pattern[i] as i64;

        if value > 0 {
            if prev & ABSORB_MASK != 0 {
                prev += value;
            }

            inputs[i] = ABSORB_MASK + value;
        }
        if value > 0 {
            if prev & ABSORB_MASK == 0 {
                prev += value;
            }
            inputs[i] = value;
        }
    }

    let mut inputs = inputs.map(|i| i as Field);
    inputs[L + 2] = domain_sep;
    inputs
}

struct SafeSponge<let L: u32> {
    state: [Field; WIDTH],
    absorb_pos: u32,
    squeeze_pos: u32,
    // expected IO pattern (compact 32-bit words)
    io_pattern: [i32; L],
    io_count: u32,
}

impl<let L: u32> SafeSponge<L> {
    pub fn start(pattern: [i32; L], domain_sep: Field) -> SafeSponge<L> {
        let mut tag = hash_to_field(tag_input(pattern, domain_sep));
        tag = tag + domain_sep;

        let mut sponge = SafeSponge::<L> {
            state: [0; WIDTH],
            absorb_pos: 0,
            squeeze_pos: 0,
            io_pattern: pattern,
            io_count: 0,
        };
        sponge.state[0] = tag;
        sponge
    }

    pub fn absorb<let S: u32>(mut self, input: [Field; S]) {
        for i in 0..S {
            if self.absorb_pos == RATE {
                //TODO: check
                self.state = poseidon_permutation(self.state);
                self.absorb_pos = 0;
            }
            let pos = self.absorb_pos + CAPACITY;
            self.state[pos] = self.state[pos] + input[i];
            self.absorb_pos += 1;
        }

        self.io_count += 1;
        self.squeeze_pos = RATE;
    }

    pub fn squeeze<let S: u32>(mut self) -> [Field; S] {
        let mut out = [0; S];
        for i in 0..S {
            if self.squeeze_pos == RATE {
                self.state = poseidon2_permutation(self.state, self.state.len());
                self.squeeze_pos = 0;
                self.absorb_pos = 0;
            }
            out[i] = self.state[self.squeeze_pos + CAPACITY];
            self.squeeze_pos += 1;
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

fn poseidon_permutation(state: [Field; WIDTH]) -> [Field; WIDTH] {
    let mut new_state = state;
    for i in 0..RATE {
        new_state[RATE] = hash_2([new_state[i], new_state[WIDTH]]);
    }
    new_state
}
