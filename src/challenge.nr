use std::hash::poseidon2_permutation;

comptime global RATE: u32 = 3;

/// Works as Poseidon Sponge to return a hash value for given transcript
struct Challenge {
    state: [Field; 4],
    cache: [Field; 3],
    cache_size: u32,
}

impl Challenge {
    // Creates a new Challenge object
    pub(crate) fn new() -> Self {
        Self { cache: [0; 3], state: [0; 4], cache_size: 0 }
    }

    pub(crate) fn perform_duplex(&mut self) {
        // add the cache into sponge state
        for i in 0..RATE {
            // We effectively zero-pad the cache by only adding to the state
            // cache that is less than the specified `cache_size`
            if i < self.cache_size {
                self.state[i] += self.cache[i];
            }
        }
        self.state = poseidon2_permutation(self.state, 4);
    }

    pub(crate) fn absorb(&mut self, input: Field) {
        if self.cache_size == RATE {
            // If we're absorbing, and the cache is full, apply the sponge permutation to compress the cache
            self.perform_duplex();
            self.cache[0] = input;
            self.cache_size = 1;
        } else {
            // If we're absorbing, and the cache is not full, add the input into the cache
            self.cache[self.cache_size] = input;
            self.cache_size += 1;
        }
    }

    pub(crate) fn squeeze(&mut self) -> Field {
        // If we're in absorb mode, apply sponge permutation to compress the cache.
        self.perform_duplex();

        // Pop one item off the top of the permutation and return it.
        self.state[0]
    }

    pub(crate) fn finish(&mut self) {
        // Overwrite state and cache to erase it
        self.state = [0; 4];
        self.cache = [0; 3];
        let res = self.squeeze();
        assert_eq(res, 0x18dfb8dc9b82229cff974efefc8df78b1ce96d9d844236b496785c698bc6732e)
    }
}

#[test]
fn test_challenge() {
    let mut challenge = Challenge::new();
    challenge.absorb(3);
    let result = challenge.squeeze();
    assert_eq(result, 0x242625e88b1aae5923a32b65ae620a8fd0830822031e994d6c241f06eb899c69);
    challenge.finish();
}
