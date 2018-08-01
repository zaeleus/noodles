#[derive(Debug)]
pub struct Reference {
    name: String,
    len: i32,
}

impl Reference {
    pub fn new(name: String, len: i32) -> Reference {
        Reference { name, len }
    }
}
