use bstr::{BStr, BString};

/// A FASTQ record definition.
///
/// A definition represents a definition line, i.e., a read name and, optionally, a description.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Definition {
    name: BString,
    description: BString,
}

impl Definition {
    /// Creates a FASTQ record definition.
    pub fn new<N, D>(name: N, description: D) -> Self
    where
        N: Into<BString>,
        D: Into<BString>,
    {
        Self {
            name: name.into(),
            description: description.into(),
        }
    }

    /// Returns the record name.
    pub fn name(&self) -> &BStr {
        self.name.as_ref()
    }

    /// Returns a mutable reference to the record name.
    pub fn name_mut(&mut self) -> &mut BString {
        &mut self.name
    }

    /// Returns the description.
    pub fn description(&self) -> &BStr {
        self.description.as_ref()
    }

    /// Returns a mutable reference to the description.
    pub fn description_mut(&mut self) -> &mut BString {
        &mut self.description
    }

    pub(crate) fn clear(&mut self) {
        self.name.clear();
        self.description.clear();
    }
}
