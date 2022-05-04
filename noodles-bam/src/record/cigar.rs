//! BAM record CIGAR.

use std::str::FromStr;

use bytes::BytesMut;
use noodles_sam as sam;
use once_cell::sync::OnceCell;

/// Lazy BAM record CIGAR.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Cigar {
    pub(crate) buf: BytesMut,
    cell: OnceCell<sam::alignment::record::Cigar>,
}

impl Cigar {
    /// Returns the number of CIGAR operations.
    pub fn len(&self) -> usize {
        self.cell
            .get()
            .map(|cigar| cigar.len())
            .unwrap_or(self.buf.len() / 4)
    }

    /// Returns whether there are any CIGAR operations.
    pub fn is_empty(&self) -> bool {
        self.cell
            .get()
            .map(|cigar| cigar.is_empty())
            .unwrap_or(self.buf.is_empty())
    }

    /// Removes all CIGAR operations.
    pub fn clear(&mut self) {}

    /// Returns an alignment record CIGAR.
    pub fn try_get(
        &self,
    ) -> Result<&sam::alignment::record::Cigar, sam::alignment::record::cigar::ParseError> {
        use crate::reader::alignment_record::get_cigar;

        self.cell.get_or_try_init(|| {
            let mut src = &self.buf[..];
            let mut cigar = sam::alignment::record::Cigar::default();
            get_cigar(&mut src, &mut cigar, self.len()).unwrap();
            Ok(cigar)
        })
    }
}

/// An error returned when a raw alignment record CIGAR fails to parse.
pub type ParseError = sam::alignment::record::cigar::ParseError;

impl FromStr for Cigar {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let sam_cigar: sam::alignment::record::Cigar = s.parse()?;
        Ok(Self::from(sam_cigar))
    }
}

impl From<sam::alignment::record::Cigar> for Cigar {
    fn from(cigar: sam::alignment::record::Cigar) -> Self {
        let cell = OnceCell::new();
        cell.set(cigar).ok();

        Self {
            buf: BytesMut::new(),
            cell,
        }
    }
}

impl TryFrom<Cigar> for sam::alignment::record::Cigar {
    type Error = ParseError;

    fn try_from(bam_cigar: Cigar) -> Result<Self, Self::Error> {
        use crate::reader::alignment_record::get_cigar;

        if let Some(sam_cigar) = bam_cigar.cell.into_inner() {
            Ok(sam_cigar)
        } else {
            let mut src = &bam_cigar.buf[..];
            let op_count = src.len() / 4;
            let mut cigar = sam::alignment::record::Cigar::default();
            get_cigar(&mut src, &mut cigar, op_count).unwrap();
            Ok(cigar)
        }
    }
}
