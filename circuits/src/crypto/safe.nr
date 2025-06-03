use crate::constants::TAG;
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
fn test_safe_sponge() {
    let pattern = [5, 4];
    let mut safe = SafeSponge::start(pattern);
    safe = safe.absorb([1, 2, 3, 4, 5]);

    let result = safe.squeeze();

    assert(result.len() == 4);
    assert(result.get(0) == 0x2ab2e2868c2490af2e13c32a0ba0852fe8d6b0979baa5831134e59029de898ab);
    assert(result.get(1) == 0x49eacd210d0bd4ae8609da50602d3144fca0ad09b77fbf3d7186068ffd374c);
    assert(result.get(2) == 0x16e39eccdecf20ef20ab765865573ffa55ad04b6e3dc44eaa2ee9a764e5f7159);
    assert(result.get(3) == 0x218adb7f25e9d014affbc444372dedcfbe9767fbf41c25412bd5bf99200be37f);

    safe.finish();
}

#[test]
fn test_safe_sponge_multiple_absorb() {
    let pattern = [3, 2];
    let mut safe = SafeSponge::start(pattern);
    safe = safe.absorb([1, 2, 3]);
    let result = safe.squeeze();

    assert(result.len() == 2);
    assert(result.get(0) == 0x0b6756778988e5a5f2127f26baf7b2056dcd9e26b7eff2716d47f3e8c4a9c38e);
    assert(result.get(1) == 0x20bd726a0eed7e98a62845891b62bee66ceab93fc41829acfef50a29e63fae0b);

    safe.finish();
}
