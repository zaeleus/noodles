/// A FASTQ record definition.
///
/// A definition represents a definition line, i.e., a read name and, optionally, a description.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Definition {
    name: Vec<u8>,
    description: Vec<u8>,
}

impl Definition {
    /// Creates a FASTQ record definition.
    pub fn new<N, D>(name: N, description: D) -> Self
    where
        N: Into<Vec<u8>>,
        D: Into<Vec<u8>>,
    {
        Self {
            name: name.into(),
            description: description.into(),
        }
    }

    /// Returns the record name.
    pub fn name(&self) -> &[u8] {
        &self.name
    }

    /// Returns a mutable reference to the record name.
    pub fn name_mut(&mut self) -> &mut Vec<u8> {
        &mut self.name
    }

    /// Returns the description.
    pub fn description(&self) -> &[u8] {
        &self.description
    }

    /// Returns a mutable reference to the description.
    pub fn description_mut(&mut self) -> &mut Vec<u8> {
        &mut self.description
    }

    pub(crate) fn clear(&mut self) {
        self.name.clear();
        self.description.clear();
    }
}
