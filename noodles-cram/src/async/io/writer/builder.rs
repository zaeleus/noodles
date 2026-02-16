use std::path::Path;

use noodles_fasta as fasta;
use tokio::{
    fs::File,
    io::{self, AsyncWrite},
};

use super::Writer;
use crate::{container::BlockContentEncoderMap, file_definition::Version, io::writer::Options};

/// An async CRAM writer builder.
#[derive(Default)]
pub struct Builder {
    reference_sequence_repository: fasta::Repository,
    options: Options,
    version_explicitly_set: bool,
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

    /// Sets the number of records per slice.
    ///
    /// The default is 10240.
    ///
    /// # Panics
    ///
    /// Panics if `records_per_slice` is 0.
    pub fn set_records_per_slice(mut self, records_per_slice: usize) -> Self {
        assert!(records_per_slice > 0, "records_per_slice must be > 0");
        self.options.records_per_slice = records_per_slice;
        self
    }

    /// Sets the number of slices per container.
    ///
    /// The default is 1.
    ///
    /// # Panics
    ///
    /// Panics if `slices_per_container` is 0.
    pub fn set_slices_per_container(mut self, slices_per_container: usize) -> Self {
        assert!(slices_per_container > 0, "slices_per_container must be > 0");
        self.options.slices_per_container = slices_per_container;
        self
    }

    /// Sets whether to embed reference sequences in slices.
    pub fn embed_reference_sequences(mut self, value: bool) -> Self {
        self.options.embed_reference_sequences = value;
        self
    }

    /// Sets whether to strip MD and NM tags from records during encoding.
    pub fn strip_md_nm(mut self, value: bool) -> Self {
        self.options.strip_md_nm = value;
        self
    }

    /// Sets the CRAM version for the output file.
    ///
    /// By default, the version is auto-detected from the block content-encoder map.
    /// Use this to explicitly request a specific version (e.g., CRAM 4.0).
    pub fn set_version(mut self, version: Version) -> Self {
        self.options.version = version;
        self.version_explicitly_set = true;
        self
    }

    /// Sets whether an external reference sequence is required.
    ///
    /// When `false`, allows writing mapped reads without a reference sequence.
    /// The default is `true`.
    pub fn set_reference_required(mut self, reference_required: bool) -> Self {
        self.options.reference_required = reference_required;
        self
    }

    /// Sets the CRAM 4.0 quality score orientation flag.
    ///
    /// When `true` (the default), quality scores are stored in alignment orientation.
    /// When `false`, quality scores are stored in original/sequencing orientation,
    /// and the reader reverses them for reverse-strand reads.
    ///
    /// This option is ignored for CRAM versions before 4.0.
    pub fn set_qs_seq_orient(mut self, qs_seq_orient: bool) -> Self {
        self.options.qs_seq_orient = qs_seq_orient;
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

        if !self.version_explicitly_set
            && uses_cram_3_1_codecs(&self.options.block_content_encoder_map)
        {
            self.options.version = Version::new(3, 1);
        }

        let records_per_container = self
            .options
            .records_per_slice
            .checked_mul(self.options.slices_per_container)
            .expect("records_per_container overflow");

        Writer {
            inner: writer,
            reference_sequence_repository: self.reference_sequence_repository,
            options: self.options,
            records: Vec::with_capacity(records_per_container),
            record_counter: 0,
        }
    }
}
