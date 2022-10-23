//! Alignment reader.

mod builder;

pub use self::builder::Builder;

use std::io::{self, Read};

use noodles_fasta as fasta;
use noodles_sam::{self as sam, alignment::Record, AlignmentReader};

/// An alignment reader.
pub struct Reader<R> {
    inner: Box<dyn AlignmentReader<R>>,
    reference_sequence_repository: fasta::Repository,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Reads and parses an alignment header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Cursor};
    /// use noodles_sam::{
    ///     self as sam,
    ///     header::record::value::{map::{self, header::Version}, Map},
    /// };
    /// use noodles_util::alignment;
    ///
    /// let data = Cursor::new(b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ");
    ///
    /// let mut reader = alignment::reader::Builder::default().build_from_reader(data)?;
    /// let header = reader.read_header()?;
    /// let actual: sam::Header = header.parse()?;
    ///
    /// let expected = sam::Header::builder()
    ///     .set_header(Map::<map::Header>::new(Version::new(1, 6)))
    ///     .build();
    ///
    /// assert_eq!(actual, expected);
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<String> {
        self.inner.read_alignment_header()
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Cursor};
    /// use noodles_sam as sam;
    /// use noodles_util::alignment;
    ///
    /// let data = Cursor::new(b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ");
    ///
    /// let mut reader = alignment::reader::Builder::default().build_from_reader(data)?;
    /// let header: sam::Header = reader.read_header()?.parse()?;
    ///
    /// let mut records = reader.records(&header);
    ///
    /// assert!(records.next().transpose()?.is_some());
    /// assert!(records.next().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn records<'a>(
        &'a mut self,
        header: &'a sam::Header,
    ) -> impl Iterator<Item = io::Result<Record>> + 'a {
        self.inner
            .alignment_records(&self.reference_sequence_repository, header)
    }
}
