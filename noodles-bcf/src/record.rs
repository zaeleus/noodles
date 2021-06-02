//! BCF record and fields.

mod convert;
pub mod value;

pub use self::value::Value;

use std::{
    convert::{TryFrom, TryInto},
    io,
    ops::{Deref, DerefMut},
};

use noodles_vcf as vcf;

/// A BCF record.
#[derive(Clone, Default, Eq, PartialEq)]
pub struct Record(Vec<u8>);

impl Record {
    pub(crate) fn resize(&mut self, new_len: usize) {
        self.0.resize(new_len, Default::default());
    }

    /// Returns the chromosome ID of the record.
    ///
    /// The chromosome ID is the index of the associated contig in the VCF header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    ///
    /// let record = bcf::Record::from(vec![
    ///     0x08, 0x00, 0x00, 0x00, // CHROM
    ///     // ...
    /// ]);
    ///
    /// assert_eq!(record.chromosome_id()?, 8);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn chromosome_id(&self) -> io::Result<i32> {
        const OFFSET: usize = 0;

        let data = &self.0[OFFSET..OFFSET + 4];

        data.try_into()
            .map(i32::from_le_bytes)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    /// Returns the start position of this record.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// use noodles_vcf as vcf;
    ///
    /// let record = bcf::Record::from(vec![
    ///     0x08, 0x00, 0x00, 0x00, // CHROM
    ///     0x0c, 0x00, 0x00, 0x00, // POS
    ///     // ...
    /// ]);
    ///
    /// assert_eq!(record.position().map(i32::from)?, 13);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn position(&self) -> io::Result<vcf::record::Position> {
        use vcf::record::Position;

        const OFFSET: usize = 4;

        let data = &self.0[OFFSET..OFFSET + 4];

        data.try_into()
            .map(i32::from_le_bytes)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|pos| {
                Position::try_from(pos + 1)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn rlen(&self) -> io::Result<i32> {
        const OFFSET: usize = 8;

        let data = &self.0[OFFSET..OFFSET + 4];

        data.try_into()
            .map(i32::from_le_bytes)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    /// Returns the end position of this record.
    ///
    /// This value is 1-based.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// use noodles_vcf as vcf;
    ///
    /// let record = bcf::Record::from(vec![
    ///     0x08, 0x00, 0x00, 0x00, // CHROM
    ///     0x0c, 0x00, 0x00, 0x00, // POS
    ///     0x05, 0x00, 0x00, 0x00, // rlen
    ///     // ...
    /// ]);
    ///
    /// assert_eq!(record.end().map(i32::from)?, 17);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn end(&self) -> io::Result<vcf::record::Position> {
        use vcf::record::Position;

        let start = self.position().map(i32::from)?;
        let len = self.rlen()?;
        let end = start + len - 1;

        Position::try_from(end).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}

impl Deref for Record {
    type Target = [u8];

    fn deref(&self) -> &[u8] {
        &self.0
    }
}

impl DerefMut for Record {
    fn deref_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl From<Vec<u8>> for Record {
    fn from(data: Vec<u8>) -> Self {
        Self(data)
    }
}
