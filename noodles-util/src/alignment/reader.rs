mod builder;

pub use self::builder::Builder;

use std::io::{self, BufReader, Read, Seek};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::{self as sam, AlignmentReader};

enum Inner<R> {
    Sam(sam::Reader<BufReader<R>>),
    Bam(bam::Reader<bgzf::Reader<R>>),
    Cram(cram::Reader<R>),
}

/// An alignment reader.
pub struct Reader<R> {
    inner: Inner<R>,
    reference_sequence_repository: fasta::Repository,
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
    ) -> impl Iterator<Item = io::Result<Box<dyn sam::AlignmentRecord>>> + 'a {
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
}
