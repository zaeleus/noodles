//! Alignment writer builder.

use std::{
    fs::File,
    io::{self, Write},
    path::Path,
};

use cram::data_container::BlockContentEncoderMap;
use noodles_bam as bam;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::Writer;
use crate::alignment::Format;

/// An alignment writer builder.
pub struct Builder {
    format: Format,
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
        self.format = format;
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
    pub fn build_from_path<P>(self, src: P) -> io::Result<Writer>
    where
        P: AsRef<Path>,
    {
        let file = File::create(src)?;
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
        let inner: Box<dyn sam::AlignmentWriter> = match self.format {
            Format::Sam => Box::new(sam::Writer::new(writer)),
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

impl Default for Builder {
    fn default() -> Self {
        Self {
            format: Format::Sam,
            reference_sequence_repository: fasta::Repository::default(),
            block_content_encoder_map: BlockContentEncoderMap::default(),
        }
    }
}
