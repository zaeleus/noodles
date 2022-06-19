mod builder;

pub use self::builder::Builder;

use std::io::{self, BufReader, Read, Seek};

use noodles_bam::{self as bam, bai};
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_cram::{self as cram, crai};
use noodles_csi as csi;
use noodles_fasta as fasta;
use noodles_sam::{self as sam, alignment::Record, AlignmentReader};

enum Inner<R> {
    Sam(sam::Reader<BufReader<R>>),
    Bam(bam::Reader<bgzf::Reader<R>>),
    Cram(cram::Reader<R>),
}

enum Index {
    Bai(bai::Index),
    Crai(crai::Index),
    Csi(csi::Index),
}

/// An alignment reader.
pub struct Reader<R> {
    inner: Inner<R>,
    reference_sequence_repository: fasta::Repository,
    index: Option<Index>,
}

impl Reader<()> {
    /// Creates an alignment reader builder.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment;
    /// let builder = alignment::Reader::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::new()
    }
}

impl<R> Reader<R>
where
    R: Read + Seek,
{
    /// Reads and parses an alignment header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Cursor};
    /// use noodles_sam::{self as sam, header::header::Version};
    /// use noodles_util::alignment;
    ///
    /// let data = Cursor::new(b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ");
    ///
    /// let mut reader = alignment::Reader::builder().build_from_reader(data)?;
    /// let actual = reader.read_header()?;
    ///
    /// let expected = sam::Header::builder()
    ///     .set_header(sam::header::header::Header::new(Version::new(1, 6)))
    ///     .build();
    ///
    /// assert_eq!(actual, expected);
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<sam::Header> {
        match &mut self.inner {
            Inner::Sam(inner) => inner.read_alignment_header(),
            Inner::Bam(inner) => inner.read_alignment_header(),
            Inner::Cram(inner) => inner.read_alignment_header(),
        }
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Cursor};
    /// use noodles_sam::{self as sam, header::header::Version};
    /// use noodles_util::alignment;
    ///
    /// let data = Cursor::new(b"@HD\tVN:1.6
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ");
    ///
    /// let mut reader = alignment::Reader::builder().build_from_reader(data)?;
    /// let header = reader.read_header()?;
    ///
    /// let mut records = reader.records(&header);
    ///
    /// assert!(records.next().transpose()?.is_some());
    /// assert!(records.next().is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records<'a>(
        &'a mut self,
        header: &'a sam::Header,
    ) -> impl Iterator<Item = io::Result<Record>> + 'a {
        match &mut self.inner {
            Inner::Sam(inner) => {
                inner.alignment_records(&self.reference_sequence_repository, header)
            }
            Inner::Bam(inner) => {
                inner.alignment_records(&self.reference_sequence_repository, header)
            }
            Inner::Cram(inner) => {
                inner.alignment_records(&self.reference_sequence_repository, header)
            }
        }
    }

    /// Returns an iterator over records that intersect the given region.
    pub fn query<'a>(
        &'a mut self,
        header: &'a sam::Header,
        region: &'a Region,
    ) -> io::Result<impl Iterator<Item = io::Result<Record>> + 'a> {
        let index = self.index.as_ref().ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "cannot query without an index")
        })?;

        let iter: Box<dyn Iterator<Item = _>> = match &mut self.inner {
            Inner::Bam(inner) => match index {
                Index::Bai(bai) => {
                    Box::new(inner.query(header.reference_sequences(), bai, region)?)
                }
                _ => todo!(),
            },
            Inner::Cram(inner) => match index {
                Index::Crai(crai) => Box::new(
                    inner
                        .query(&self.reference_sequence_repository, header, crai, region)?
                        .map(|result| {
                            result.and_then(|record| record.try_into_alignment_record(header))
                        }),
                ),
                _ => todo!(),
            },
            _ => todo!(),
        };

        Ok(iter)
    }
}
