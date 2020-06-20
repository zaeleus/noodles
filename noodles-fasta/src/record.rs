pub mod definition;

pub use self::definition::Definition;

pub struct Record {
    definition: Definition,
    sequence: Vec<u8>,
}

impl Record {
    pub fn new(definition: Definition, sequence: Vec<u8>) -> Self {
        Self {
            definition,
            sequence,
        }
    }

    pub fn name(&self) -> &str {
        &self.definition.reference_sequence_name()
    }

    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }
}
