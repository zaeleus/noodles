use std::io;

/// A raw SAM record sequence.
#[derive(Debug, Eq, PartialEq)]
pub struct Sequence<'a>(&'a [u8]);

impl<'a> Sequence<'a> {
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self(buf)
    }

    /// Returns whether the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of bases in the sequence.
    pub fn len(&self) -> usize {
        self.0.len()
    }
}

impl<'a> AsRef<[u8]> for Sequence<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> TryInto<crate::record::Sequence> for Sequence<'a> {
    type Error = io::Error;

    fn try_into(self) -> Result<crate::record::Sequence, Self::Error> {
        use crate::reader::record::parse_sequence;

        let mut sequence = crate::record::Sequence::default();

        if !self.is_empty() {
            parse_sequence(self.as_ref(), &mut sequence)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        }

        Ok(sequence)
    }
}
