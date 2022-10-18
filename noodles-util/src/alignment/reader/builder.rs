use std::{
    fs::File,
    io::{self, BufRead, BufReader, Read, Seek},
    path::Path,
};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;
use sam::AlignmentReader;

use super::Reader;
use crate::alignment::Format;

/// An alignment reader builder.
#[derive(Default)]
pub struct Builder {
    format: Option<Format>,
    reference_sequence_repository: fasta::Repository,
}

impl Builder {
    /// Sets the format of the input.
    ///
    /// By default, the format is autodetected on [`build`]. This can be used to override it.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::alignment::{self, Format};
    /// let builder = alignment::reader::Builder::default().set_format(Format::Sam);
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
    /// use noodles_fasta as fasta;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let repository = fasta::Repository::default();
    ///
    /// let builder = alignment::reader::Builder::default()
    ///     .set_reference_sequence_repository(repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Builds an alignment reader from a path.
    ///
    /// By default, the format will be autodetected. This can be overridden by using
    /// [`set_format`]. An associated index will also attempt to be loaded.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_util::alignment;
    /// let reader = alignment::reader::Builder::default().build_from_path("sample.bam")?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, path: P) -> io::Result<Reader<Box<dyn BufRead>>>
    where
        P: AsRef<Path>,
    {
        let file = File::open(path)?;
        self.build_from_reader(file)
    }

    /// Builds an alignment reader from a reader.
    ///
    /// By default, the format will be autodetected. This can be overridden by using
    /// [`set_format`].
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment;
    /// let reader = alignment::reader::Builder::default().build_from_reader(io::empty())?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_reader<R>(self, mut reader: R) -> io::Result<Reader<Box<dyn BufRead>>>
    where
        R: Read + Seek + 'static,
    {
        let format = self
            .format
            .map(Ok)
            .unwrap_or_else(|| detect_format(&mut reader))?;

        let inner: Box<dyn AlignmentReader<_>> = match format {
            Format::Sam => {
                let inner: Box<dyn BufRead> = Box::new(BufReader::new(reader));
                Box::new(sam::Reader::from(inner))
            }
            Format::Bam => {
                let inner: Box<dyn BufRead> = Box::new(bgzf::Reader::new(reader));
                Box::new(bam::Reader::from(inner))
            }
            Format::Cram => {
                let inner: Box<dyn BufRead> = Box::new(BufReader::new(reader));
                Box::new(cram::Reader::new(inner))
            }
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
