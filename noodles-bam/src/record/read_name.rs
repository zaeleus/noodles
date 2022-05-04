//! BAM record read name.

use std::str::FromStr;

use bytes::BytesMut;
use noodles_sam as sam;
use once_cell::sync::OnceCell;

/// A lazy BAM record read name.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ReadName {
    pub(crate) buf: BytesMut,
    cell: OnceCell<sam::alignment::record::ReadName>,
}

impl ReadName {
    /// Returns an alignment record read name.
    pub fn try_get(
        &self,
    ) -> Result<&sam::alignment::record::ReadName, sam::alignment::record::read_name::ParseError>
    {
        self.cell
            .get_or_try_init(|| sam::alignment::record::ReadName::try_new(&self.buf[..]))
    }
}

/// An error returned when a raw alignment record read name fails to parse.
pub type ParseError = sam::alignment::record::read_name::ParseError;

impl FromStr for ReadName {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let sam_read_name: sam::alignment::record::ReadName = s.parse()?;
        Ok(Self::from(sam_read_name))
    }
}

impl From<sam::alignment::record::ReadName> for ReadName {
    fn from(read_name: sam::alignment::record::ReadName) -> Self {
        let cell = OnceCell::new();
        cell.set(read_name).ok();

        Self {
            buf: BytesMut::new(),
            cell,
        }
    }
}

impl From<ReadName> for sam::alignment::record::ReadName {
    fn from(bam_read_name: ReadName) -> Self {
        if let Some(sam_read_name) = bam_read_name.cell.into_inner() {
            sam_read_name
        } else {
            let buf = bam_read_name.buf;
            Self::try_new(&buf[..]).unwrap()
        }
    }
}
