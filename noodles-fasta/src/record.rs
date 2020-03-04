#[derive(Clone, Debug, Default)]
pub struct Record {
    name: String,
    sequence: Vec<u8>,
}

impl Record {
    pub fn new(name: String, sequence: Vec<u8>) -> Self {
        Self { name, sequence }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub(crate) fn name_mut(&mut self) -> &mut String {
        &mut self.name
    }

    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    pub(crate) fn sequence_mut(&mut self) -> &mut Vec<u8> {
        &mut self.sequence
    }

    pub(crate) fn clear(&mut self) {
        self.name.clear();
        self.sequence.clear();
    }
}
