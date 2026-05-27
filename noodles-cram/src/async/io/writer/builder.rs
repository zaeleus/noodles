use std::path::Path;

use noodles_fasta as fasta;
use tokio::{
    fs::File,
    io::{self, AsyncWrite},
};

use super::Writer;
use crate::{
    container::BlockContentEncoderMap,
    file_definition::Version,
    io::writer::{Context, DEFAULT_SLICES_PER_CONTAINER, builder::Preset},
};

/// An async CRAM writer builder.
#[derive(Default)]
pub struct Builder {
    context: Context,
    preset: Preset,
    block_content_encoder_map: Option<BlockContentEncoderMap>,
}

impl Builder {
    /// Sets the reference sequence repository.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::r#async::io::writer::Builder;
    /// use noodles_fasta as fasta;
    ///
    /// let repository = fasta::Repository::default();
    /// let builder = Builder::default()
    ///     .set_reference_sequence_repository(repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.context.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Sets whether to preserve read names.
    ///
    /// If `false`, read names are discarded.
    ///
    /// The default is `true`.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::r#async::io::writer::Builder;
    /// let builder = Builder::default().preserve_read_names(false);
    /// ```
    pub fn preserve_read_names(mut self, value: bool) -> Self {
        self.context.preserve_read_names = value;
        self
    }

    /// Sets whether to encode alignment start positions as deltas.
    ///
    /// If `false`, record alignment start positions are written with their actual values.
    ///
    /// The default is `true`.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::r#async::io::writer::Builder;
    /// let builder = Builder::default()
    ///     .encode_alignment_start_positions_as_deltas(false);
    /// ```
    pub fn encode_alignment_start_positions_as_deltas(mut self, value: bool) -> Self {
        self.context.encode_alignment_start_positions_as_deltas = value;
        self
    }

    /// Sets the block content-encoder map.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::{r#async::io::writer::Builder, container::BlockContentEncoderMap};
    ///
    /// let block_content_encoder_map = BlockContentEncoderMap::default();
    /// let builder = Builder::default()
    ///     .set_block_content_encoder_map(block_content_encoder_map);
    /// ```
    pub fn set_block_content_encoder_map(
        mut self,
        block_content_encoder_map: BlockContentEncoderMap,
    ) -> Self {
        self.block_content_encoder_map = Some(block_content_encoder_map);
        self
    }

    /// Builds an async CRAM writer from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_cram::r#async::io::writer::Builder;
    /// let writer = Builder::default().build_from_path("out.cram").await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn build_from_path<P>(self, dst: P) -> io::Result<Writer<File>>
    where
        P: AsRef<Path>,
    {
        File::create(dst)
            .await
            .map(|file| self.build_from_writer(file))
    }

    /// Builds an async CRAM writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// use tokio::io;
    /// let writer = cram::r#async::io::writer::Builder::default().build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<W>(mut self, writer: W) -> Writer<W>
    where
        W: AsyncWrite + Unpin,
    {
        use crate::io::writer::builder::uses_cram_3_1_codecs;

        let block_content_encoder_map = self
            .block_content_encoder_map
            .unwrap_or_else(|| self.preset.block_content_encoder_map());

        if uses_cram_3_1_codecs(&self.context.block_content_encoder_map) {
            self.context.version = Version::new(3, 1);
        }

        let records_per_slice = self.preset.records_per_slice();
        let records_per_container = DEFAULT_SLICES_PER_CONTAINER * records_per_slice;

        self.context.records_per_slice = records_per_slice;
        self.context.block_content_encoder_map = block_content_encoder_map;

        Writer {
            inner: writer,
            context: self.context,
            records: Vec::with_capacity(records_per_container),
            record_counter: 0,
        }
    }
}
