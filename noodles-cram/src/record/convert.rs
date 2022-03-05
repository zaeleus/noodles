use std::{io, str};

use noodles_bam::record::ReferenceSequenceId;
use noodles_sam::{
    self as sam,
    record::{Data, QualityScores, Sequence},
    AlignmentRecord,
};

use super::{resolve::resolve_features, Record, Tag};

impl Record {
    /// Converts this CRAM record to a SAM record.
    ///
    /// This assumes this record is fully resolved.
    pub fn try_into_sam_record(
        &self,
        reference_sequences: &sam::header::ReferenceSequences,
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

        if let Some(alignment_start) = self.alignment_start() {
            builder = builder.set_position(alignment_start);
        }

        if let Some(mapping_quality) = self.mapping_quality() {
            builder = builder.set_mapping_quality(mapping_quality);
        }

        if !self.bam_flags().is_unmapped() {
            let cigar = resolve_features(self.features(), self.read_length() as i32);
            builder = builder.set_cigar(cigar);
        }

        if let Some(mate_reference_sequence_name) = get_reference_sequence_name(
            reference_sequences,
            self.next_fragment_reference_sequence_id(),
        )? {
            builder = builder.set_mate_reference_sequence_name(mate_reference_sequence_name);
        }

        if let Some(next_mate_alignment_start) = self.next_mate_alignment_start() {
            builder = builder.set_mate_position(next_mate_alignment_start);
        }

        builder = builder.set_template_length(self.template_size());

        if self.read_length() > 0 {
            let sequence = bytes_to_sequence(self.bases())?;
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
        .map(usize::from)
        .map(|id| {
            reference_sequences
                .get_index(id)
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
        let sam_tag = tag
            .key()
            .tag()
            .try_into()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let sam_value = tag.value().clone().into();

        let field = Field::new(sam_tag, sam_value);
        fields.push(field);
    }

    Data::try_from(fields).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
}
