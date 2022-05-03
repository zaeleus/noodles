//! BAM record data.

use std::{io, str::FromStr};

use bytes::BytesMut;
use noodles_sam as sam;
use once_cell::sync::OnceCell;

/// Lazy BAM record data.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Data {
    pub(crate) buf: BytesMut,
    cell: OnceCell<sam::alignment::record::Data>,
}

impl Data {
    /// Returns whether there are any data fields.
    pub fn is_empty(&self) -> bool {
        self.cell
            .get()
            .map(|data| data.is_empty())
            .unwrap_or(self.buf.is_empty())
    }

    /// Removes all data fields.
    pub fn clear(&mut self) {
        self.buf.clear();
        self.cell.take();
    }

    /// Returns alignment record data.
    pub fn try_get(&self) -> io::Result<&sam::alignment::record::Data> {
        use crate::reader::alignment_record::get_data;

        self.cell.get_or_try_init(|| {
            let mut src = &self.buf[..];
            let mut data = sam::alignment::record::Data::default();
            get_data(&mut src, &mut data)?;
            Ok(data)
        })
    }
}

/// An error returned when raw alignment record data fails to parse.
pub type ParseError = sam::alignment::record::data::ParseError;

impl FromStr for Data {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let sam_data: sam::alignment::record::Data = s.parse()?;
        Ok(Self::from(sam_data))
    }
}

impl From<sam::alignment::record::Data> for Data {
    fn from(sam_data: sam::alignment::record::Data) -> Self {
        let cell = OnceCell::new();
        cell.set(sam_data).ok();

        Self {
            buf: BytesMut::new(),
            cell,
        }
    }
}

impl From<Data> for sam::alignment::record::Data {
    fn from(bam_data: Data) -> Self {
        use crate::reader::alignment_record::get_data;

        if let Some(sam_data) = bam_data.cell.into_inner() {
            sam_data
        } else {
            let mut src = &bam_data.buf[..];
            let mut data = sam::alignment::record::Data::default();
            get_data(&mut src, &mut data).unwrap();
            data
        }
    }
}
