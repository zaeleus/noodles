use noodles_fasta as fasta;
use tokio::io::AsyncWrite;

use super::Writer;
use crate::{
    data_container::BlockContentEncoderMap, file_definition::Version, io::writer::Options,
    DataContainer,
};

/// An async CRAM writer builder.
#[derive(Default)]
pub struct Builder {
    reference_sequence_repository: fasta::Repository,
    options: Options,
}

impl Builder {
    /// Sets the reference sequence repository.
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Sets whether to preserve read names.
    ///
    /// If `false`, read names are discarded.
    ///
    /// The default is `true`.
    pub fn preserve_read_names(mut self, value: bool) -> Self {
        self.options.preserve_read_names = value;
        self
    }

    /// Sets whether to encode alignment start positions as deltas.
    ///
    /// If `false`, record alignment start positions are written with their actual values.
    ///
    /// The default is `true`.
    pub fn encode_alignment_start_positions_as_deltas(mut self, value: bool) -> Self {
        self.options.encode_alignment_start_positions_as_deltas = value;
        self
    }

    /// Sets the block content-encoder map.
    pub fn set_block_content_encoder_map(mut self, map: BlockContentEncoderMap) -> Self {
        self.options.block_content_encoder_map = map;
        self
    }

    /// Builds an async CRAM writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// use tokio::io;
    /// let writer = cram::r#async::io::writer::Builder::default().build_with_writer(io::sink());
    /// ```
    pub fn build_with_writer<W>(mut self, writer: W) -> Writer<W>
    where
        W: AsyncWrite + Unpin,
    {
        use crate::io::writer::builder::uses_cram_3_1_codecs;

        if uses_cram_3_1_codecs(&self.options.block_content_encoder_map) {
            self.options.version = Version::new(3, 1);
        }

        Writer {
            inner: writer,
            reference_sequence_repository: self.reference_sequence_repository,
            options: self.options,
            data_container_builder: DataContainer::builder(0),
            record_counter: 0,
        }
    }
}
