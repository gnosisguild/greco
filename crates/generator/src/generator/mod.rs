/// @todo implement generator here
pub struct Generator {
    x: u64
}

impl Generator {
    pub fn new() -> Self {
       Self { x: 0 }
    }

    pub fn generate(&mut self) {
        self.x += 1;
    }
}