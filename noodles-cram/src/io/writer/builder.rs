use std::io::Write;

use noodles_fasta as fasta;

use super::{Options, Writer};
use crate::{
    codecs::Encoder, data_container::BlockContentEncoderMap, file_definition::Version,
    DataContainer,
};

/// A CRAM writer builder.
#[derive(Default)]
pub struct Builder {
    reference_sequence_repository: fasta::Repository,
    options: Options,
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

    /// Sets the block content-encoder map.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::{data_container::BlockContentEncoderMap, io::writer::Builder};
    ///
    /// let block_content_encoder_map = BlockContentEncoderMap::default();
    /// let builder = Builder::default()
    ///     .set_block_content_encoder_map(block_content_encoder_map);
    /// ```
    pub fn set_block_content_encoder_map(mut self, map: BlockContentEncoderMap) -> Self {
        self.options.block_content_encoder_map = map;
        self
    }

    /// Builds a CRAM writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::io::writer::Builder;
    /// let writer = Builder::default().build_with_writer(Vec::new());
    /// ```
    pub fn build_with_writer<W>(mut self, writer: W) -> Writer<W>
    where
        W: Write,
    {
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

pub fn uses_cram_3_1_codecs(block_content_encoder_map: &BlockContentEncoderMap) -> bool {
    fn is_cram_3_1_codec(encoder: &Encoder) -> bool {
        matches!(
            encoder,
            Encoder::RansNx16(_) | Encoder::AdaptiveArithmeticCoding(_) | Encoder::NameTokenizer
        )
    }

    if let Some(encoder) = block_content_encoder_map.core_data_encoder() {
        if is_cram_3_1_codec(encoder) {
            return true;
        }
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
