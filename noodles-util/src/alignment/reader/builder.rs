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
    reference_sequence_repository: Option<fasta::Repository>,
}

impl<R> Builder<R>
where
    R: Read + Seek + 'static,
{
    pub(super) fn new(inner: R) -> Self {
        Self {
            inner,
            reference_sequence_repository: None,
            format: None,
        }
    }

    /// Sets the format of the input.
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Sets the reference sequence repository.
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = Some(reference_sequence_repository);
        self
    }

    /// Builds an alignment reader.
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

        let reference_sequence_repository = self.reference_sequence_repository.unwrap_or_default();

        Ok(Reader {
            inner,
            reference_sequence_repository,
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
