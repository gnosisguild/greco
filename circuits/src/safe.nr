use poseidon::poseidon::bn254::hash_2;
use std::hash::hash_to_field;
use std::hash::poseidon2_permutation;

global RATE: u32 = 4;
global CAPACITY: u32 = 1;
global WIDTH: u32 = 5;
global ABSORB_MASK: u32 = 0x80000_0000;

pub enum Pattern {
    Absorb(u32),
    Squeeze(u32),
}

impl Pattern {
    pub fn length(self) -> u32 {
        match self {
            Pattern::Absorb(len) => len,
            Pattern::Squeeze(len) => len,
        }
    }
}

fn tag_input<let L: u32>(pattern: [Pattern; L], domain_sep: u32) -> Vec<u32>{
    assert(pattern[0] == Pattern::Absorb(_));
    assert(pattern[L] == Pattern::Squeeze(_));
    let mut input = Vec::new();
    input.push(ABSORB_MASK);

    for i in 0..pattern.len(){
        let l = input.len();
        let mut prev = input.get(l - 1);

        match pattern {
            Pattern::Absorb(len) => if prev & ABSORB_MASK != 0 { prev += len; },
            Pattern::Absorb(len) => input.push(ABSORB_MASK + len),
            Pattern::Squeeze(len) => if prev & ABSORB_MASK == 0 { prev += len; },
            Pattern::Squeeze(len) => input.push(len),
        }
    }

    input.push(domain_sep);
}


struct SafeSponge<let L: u32> {
    state: [Field; WIDTH],
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
