#[derive(Clone, Debug)]
pub struct Reference {
    name: String,
    len: i32,
}

impl Reference {
    pub fn new(name: String, len: i32) -> Reference {
        Reference { name, len }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn len(&self) -> i32 {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}
