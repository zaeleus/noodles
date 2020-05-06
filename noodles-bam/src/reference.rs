#[derive(Clone, Debug)]
pub struct Reference {
    name: String,
    len: u32,
}

impl Reference {
    pub fn new(name: String, len: u32) -> Self {
        Self { name, len }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn len(&self) -> u32 {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}
