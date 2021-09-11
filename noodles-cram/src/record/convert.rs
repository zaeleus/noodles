use std::{convert::TryFrom, io, str};

use noodles_bam::record::ReferenceSequenceId;
use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    record::{Data, Position, QualityScores, Sequence},
};

use crate::data_container::CompressionHeader;

use super::{
    resolve::{resolve_bases, resolve_features},
    Record, Tag,
};

impl Record {
    /// Converts this CRAM record to a SAM record.
    pub fn try_into_sam_record(
        &self,
        reference_assembly: &[fasta::Record],
        reference_sequences: &sam::header::ReferenceSequences,
        compression_header: &CompressionHeader,
    ) -> io::Result<sam::Record> {
        let mut builder = sam::Record::builder();

        let raw_read_name = str::from_utf8(self.read_name())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        if raw_read_name != "*" {
            let read_name = raw_read_name
                .parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

            builder = builder.set_read_name(read_name);
        }

        builder = builder.set_flags(self.bam_flags());

        if let Some(reference_sequence_name) =
            get_reference_sequence_name(reference_sequences, self.reference_sequence_id())?
        {
            builder = builder.set_reference_sequence_name(reference_sequence_name);
        }

        if self.alignment_start() > 0 {
            let position = Position::try_from(self.alignment_start())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

            builder = builder.set_position(position);
        }

        builder = builder.set_mapping_quality(self.mapping_quality());

        let cigar = resolve_features(self.features(), self.read_length() as i32);
        builder = builder.set_cigar(cigar);

        if let Some(mate_reference_sequence_name) = get_reference_sequence_name(
            reference_sequences,
            self.next_fragment_reference_sequence_id(),
        )? {
            builder = builder.set_mate_reference_sequence_name(mate_reference_sequence_name);
        }

        if self.next_mate_alignment_start() > 0 {
            let mate_position = Position::try_from(self.next_mate_alignment_start())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

            builder = builder.set_mate_position(mate_position);
        }

        builder = builder.set_template_length(self.template_size());

        if self.reference_sequence_id().is_some() && self.read_length() > 0 {
            let reference_sequence_record = self
                .reference_sequence_id()
                .map(i32::from)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "missing reference sequence ID")
                })
                .and_then(|id| {
                    usize::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                })
                .and_then(|id| {
                    reference_assembly.get(id).ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "missing reference sequence")
                    })
                })?;

            let raw_bases = resolve_bases(
                reference_sequence_record,
                compression_header,
                self.features(),
                self.alignment_start(),
                self.read_length() as usize,
            );

            let sequence = bytes_to_sequence(&raw_bases)?;
            builder = builder.set_sequence(sequence);
        }

        if !self.quality_scores().is_empty() {
            let quality_scores = bytes_to_quality_scores(self.quality_scores())?;
            builder = builder.set_quality_scores(quality_scores);
        }

        if !self.tags().is_empty() {
            let data = tags_to_data(self.tags())?;
            builder = builder.set_data(data);
        }

        builder
            .build()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}

fn get_reference_sequence_name(
    reference_sequences: &sam::header::ReferenceSequences,
    reference_sequence_id: Option<ReferenceSequenceId>,
) -> io::Result<Option<sam::record::ReferenceSequenceName>> {
    reference_sequence_id
        .map(i32::from)
        .map(|id| {
            reference_sequences
                .get_index(id as usize)
                .and_then(|(_, rs)| rs.name().parse().ok())
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "invalid reference sequence ID")
                })
        })
        .transpose()
}

fn bytes_to_sequence(data: &[u8]) -> io::Result<Sequence> {
    use sam::record::sequence::Base;

    data.iter()
        .copied()
        .map(|b| Base::try_from(b as char))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        .map(Sequence::from)
}

fn bytes_to_quality_scores(data: &[u8]) -> io::Result<QualityScores> {
    use sam::record::quality_scores::Score;

    data.iter()
        .copied()
        .map(Score::try_from)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        .map(QualityScores::from)
}

fn tags_to_data(tags: &[Tag]) -> io::Result<Data> {
    use sam::record::data::Field;

    let mut fields = Vec::with_capacity(tags.len());

    for tag in tags {
        let raw_tag = tag.key().tag();
        let sam_tag = str::from_utf8(&raw_tag)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            })?;

        let bam_value = tag.value();
        let sam_value = bam_value.clone().into();

        let field = Field::new(sam_tag, sam_value);
        fields.push(field);
    }

    Data::try_from(fields).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
}
