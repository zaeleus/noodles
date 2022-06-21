use std::io;

use noodles_sam as sam;

/// Raw BAM record data.
#[derive(Debug, Eq, PartialEq)]
pub struct Data<'a>(&'a [u8]);

impl<'a> Data<'a> {
    pub(super) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }

    /// Returns whether there are any fields.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl<'a> AsRef<[u8]> for Data<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> TryFrom<Data<'a>> for sam::record::Data {
    type Error = io::Error;

    fn try_from(bam_data: Data<'a>) -> Result<Self, Self::Error> {
        use crate::reader::record::get_data;

        let mut src = bam_data.0;
        let mut sam_data = Self::default();
        get_data(&mut src, &mut sam_data)?;

        Ok(sam_data)
    }
}
