//! Alignment writer builder.

use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use cram::data_container::BlockContentEncoderMap;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::Writer;
use crate::alignment::Format;

/// An alignment writer builder.
#[derive(Default)]
pub struct Builder {
    format: Option<Format>,
    reference_sequence_repository: fasta::Repository,
    block_content_encoder_map: BlockContentEncoderMap,
}

impl Builder {
    /// Sets the format of the output.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::alignment::{self, Format};
    /// let builder = alignment::writer::Builder::default().set_format(Format::Sam);
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
    /// let builder = alignment::writer::Builder::default()
    ///     .set_reference_sequence_repository(repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Sets the block content-encoder map.
    ///
    /// This is only used when the output format is CRAM.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram::data_container::BlockContentEncoderMap;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let builder = alignment::writer::Builder::default()
    ///     .set_format(Format::Cram)
    ///     .set_block_content_encoder_map(BlockContentEncoderMap::default());
    /// ```
    pub fn set_block_content_encoder_map(
        mut self,
        block_content_encoder_map: BlockContentEncoderMap,
    ) -> Self {
        self.block_content_encoder_map = block_content_encoder_map;
        self
    }

    /// Builds an alignment writer from a path.
    ///
    /// If the format is not set, it is detected from the path extension.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let writer = alignment::writer::Builder::default()
    ///     .set_format(Format::Sam)
    ///     .build_from_path("out.sam")?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_path<P>(mut self, src: P) -> io::Result<Writer>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        if self.format.is_none() {
            self.format = detect_format_from_path_extension(src);
        }

        let file = File::create(src).map(BufWriter::new)?;
        Ok(self.build_from_writer(file))
    }

    /// Builds an alignment writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let writer = alignment::writer::Builder::default()
    ///     .set_format(Format::Sam)
    ///     .build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> Writer
    where
        W: Write + 'static,
    {
        let format = self.format.unwrap_or(Format::Sam);

        let inner: Box<dyn sam::AlignmentWriter> = match format {
            Format::Sam => Box::new(sam::Writer::new(writer)),
            Format::SamGz => Box::new(sam::Writer::new(bgzf::Writer::new(writer))),
            Format::Bam => Box::new(bam::Writer::new(writer)),
            Format::Cram => Box::new(
                cram::writer::Builder::default()
                    .set_reference_sequence_repository(self.reference_sequence_repository)
                    .set_block_content_encoder_map(self.block_content_encoder_map)
                    .build_with_writer(writer),
            ),
        };

        Writer { inner }
    }
}

fn detect_format_from_path_extension<P>(path: P) -> Option<Format>
where
    P: AsRef<Path>,
{
    let ext = path.as_ref().file_name().and_then(|ext| ext.to_str())?;

    if ext.ends_with("sam") {
        Some(Format::Sam)
    } else if ext.ends_with("sam.gz") || ext.ends_with("sam.bgz") {
        Some(Format::SamGz)
    } else if ext.ends_with("bam") {
        Some(Format::Bam)
    } else if ext.ends_with("cram") {
        Some(Format::Cram)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_format_from_path_extension() {
        assert_eq!(
            detect_format_from_path_extension("out.sam"),
            Some(Format::Sam)
        );
        assert_eq!(
            detect_format_from_path_extension("out.sam.gz"),
            Some(Format::SamGz)
        );
        assert_eq!(
            detect_format_from_path_extension("out.bam"),
            Some(Format::Bam)
        );
        assert_eq!(
            detect_format_from_path_extension("out.cram"),
            Some(Format::Cram)
        );

        assert!(detect_format_from_path_extension("out.fa").is_none());
    }
}
