pub type Id = [u8; 2];

#[derive(Clone, Debug)]
pub struct Tag {
    id: Id,
    data: Vec<u8>,
}

impl Tag {
    pub fn new(id: Id, data: Vec<u8>) -> Self {
        Self { id, data }
    }
}
