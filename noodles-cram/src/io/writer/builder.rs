use std::{
    fs::File,
    io::{self, Write},
    path::Path,
};

use noodles_fasta as fasta;

use super::{Options, Writer};
use crate::{codecs::Encoder, container::BlockContentEncoderMap, file_definition::Version};

/// A CRAM writer builder.
#[derive(Default)]
pub struct Builder {
    reference_sequence_repository: fasta::Repository,
    options: Options,
    version_explicitly_set: bool,
}

impl Builder {
    /// Sets the reference sequence repository.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::io::writer::Builder;
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
        self.reference_sequence_repository = reference_sequence_repository;
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
    /// use noodles_cram::io::writer::Builder;
    /// let builder = Builder::default().preserve_read_names(false);
    /// ```
    pub fn preserve_read_names(mut self, value: bool) -> Self {
        self.options.preserve_read_names = value;
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
    /// use noodles_cram::io::writer::Builder;
    /// let builder = Builder::default()
    ///     .encode_alignment_start_positions_as_deltas(false);
    /// ```
    pub fn encode_alignment_start_positions_as_deltas(mut self, value: bool) -> Self {
        self.options.encode_alignment_start_positions_as_deltas = value;
        self
    }

    /// Sets the CRAM version for the output file.
    ///
    /// By default, the version is auto-detected from the block content-encoder map.
    /// Use this to explicitly request a specific version (e.g., CRAM 4.0).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::{file_definition::Version, io::writer::Builder};
    /// let builder = Builder::default().set_version(Version::new(4, 0));
    /// ```
    pub fn set_version(mut self, version: Version) -> Self {
        self.options.version = version;
        self.version_explicitly_set = true;
        self
    }

    /// Sets the block content-encoder map.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::{container::BlockContentEncoderMap, io::writer::Builder};
    ///
    /// let block_content_encoder_map = BlockContentEncoderMap::default();
    /// let builder = Builder::default()
    ///     .set_block_content_encoder_map(block_content_encoder_map);
    /// ```
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
    ///
    /// When enabled, reference subsequences are stored directly in the CRAM file,
    /// removing the need for an external reference file during decoding.
    ///
    /// The default is `false`.
    pub fn embed_reference_sequences(mut self, value: bool) -> Self {
        self.options.embed_reference_sequences = value;
        self
    }

    /// Sets whether to strip MD and NM tags from records during encoding.
    ///
    /// When enabled, MD and NM tags are omitted from the output since they can
    /// be reconstructed from CRAM features and the reference.
    ///
    /// The default is `false`.
    pub fn strip_md_nm(mut self, value: bool) -> Self {
        self.options.strip_md_nm = value;
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

    /// Builds a CRAM writer from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_cram::io::writer::Builder;
    /// let writer = Builder::default().build_from_path("out.cram")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, dst: P) -> io::Result<Writer<File>>
    where
        P: AsRef<Path>,
    {
        File::create(dst).map(|file| self.build_from_writer(file))
    }

    /// Builds a CRAM writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::io::writer::Builder;
    /// let writer = Builder::default().build_from_writer(Vec::new());
    /// ```
    pub fn build_from_writer<W>(mut self, writer: W) -> Writer<W>
    where
        W: Write,
    {
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

pub fn uses_cram_3_1_codecs(block_content_encoder_map: &BlockContentEncoderMap) -> bool {
    fn is_cram_3_1_codec(encoder: &Encoder) -> bool {
        matches!(
            encoder,
            Encoder::RansNx16(_) | Encoder::AdaptiveArithmeticCoding(_) | Encoder::NameTokenizer
        )
    }

    if let Some(encoder) = block_content_encoder_map.core_data_encoder()
        && is_cram_3_1_codec(encoder)
    {
        return true;
    }

    block_content_encoder_map
        .data_series_encoders()
        .iter()
        .chain(block_content_encoder_map.tag_values_encoders().values())
        .flatten()
        .any(is_cram_3_1_codec)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uses_cram_3_1_codecs() {
        use crate::codecs::rans_nx16::Flags;

        let block_content_encoder_map = BlockContentEncoderMap::default();
        assert!(!uses_cram_3_1_codecs(&block_content_encoder_map));

        let block_content_encoder_map = BlockContentEncoderMap::builder()
            .set_core_data_encoder(Some(Encoder::RansNx16(Flags::empty())))
            .build();
        assert!(uses_cram_3_1_codecs(&block_content_encoder_map));
    }
}
