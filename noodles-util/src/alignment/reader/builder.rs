use std::io::{self, BufReader, Read, Seek};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::Reader;
use crate::alignment::Format;

/// An alignment reader builder.
pub struct Builder<R> {
    inner: R,
    format: Option<Format>,
    reference_sequence_repository: fasta::Repository,
}

impl<R> Builder<R>
where
    R: Read + Seek + 'static,
{
    pub(super) fn new(inner: R) -> Self {
        Self {
            inner,
            reference_sequence_repository: fasta::Repository::default(),
            format: None,
        }
    }

    /// Sets the format of the input.
    ///
    /// By default, the format is autodetected on [`build`]. This can be used to override it.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment::{self, Format};
    /// let builder = alignment::Reader::builder(io::empty()).set_format(Format::Sam);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Sets the reference sequence repository.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta as fasta;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let repository = fasta::Repository::default();
    ///
    /// let builder = alignment::Reader::builder(io::empty())
    ///     .set_reference_sequence_repository(repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Builds an alignment reader.
    ///
    /// By default, the format will be autodetected. This can be overridden by using
    /// [`set_format`].
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment;
    /// let reader = alignment::Reader::builder(io::empty()).build()?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build(mut self) -> io::Result<Reader> {
        let format = self
            .format
            .map(Ok)
            .unwrap_or_else(|| detect_format(&mut self.inner))?;

        let inner: Box<dyn sam::AlignmentReader> = match format {
            Format::Sam => Box::new(sam::Reader::new(BufReader::new(self.inner))),
            Format::Bam => Box::new(bam::Reader::new(self.inner)),
            Format::Cram => Box::new(cram::Reader::new(self.inner)),
        };

        Ok(Reader {
            inner,
            reference_sequence_repository: self.reference_sequence_repository,
        })
    }
}

fn detect_format<R>(reader: &mut R) -> io::Result<Format>
where
    R: Read + Seek,
{
    const CRAM_MAGIC_NUMBER: [u8; 4] = [b'C', b'R', b'A', b'M'];
    const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];
    const BAM_MAGIC_NUMBER: [u8; 4] = [b'B', b'A', b'M', 0x01];

    let mut buf = [0; 4];
    reader.read_exact(&mut buf).ok();
    reader.rewind()?;

    if buf == CRAM_MAGIC_NUMBER {
        return Ok(Format::Cram);
    }

    if buf[..2] == GZIP_MAGIC_NUMBER {
        let mut reader = bgzf::Reader::new(reader);
        reader.read_exact(&mut buf).ok();
        reader.get_mut().rewind()?;

        if buf == BAM_MAGIC_NUMBER {
            return Ok(Format::Bam);
        }
    }

    Ok(Format::Sam)
}
