use std::io;

use noodles_bam::record::ReferenceSequenceId;
use noodles_sam::{self as sam, record::Data, AlignmentRecord};

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

        if let Some(read_name) = self.read_name() {
            builder = builder.set_read_name(read_name.clone());
        }

        builder = builder.set_flags(self.bam_flags());

        if let Some(reference_sequence_name) =
            get_reference_sequence_name(reference_sequences, self.reference_sequence_id())?
        {
            builder = builder.set_reference_sequence_name(reference_sequence_name);
        }

        if let Some(alignment_start) = self.alignment_start() {
            let alignment_start =
                sam::record::Position::try_from(usize::from(alignment_start) as i32).unwrap();
            builder = builder.set_position(alignment_start);
        }

        if let Some(mapping_quality) = self.mapping_quality() {
            builder = builder.set_mapping_quality(mapping_quality);
        }

        if !self.bam_flags().is_unmapped() {
            let cigar = resolve_features(self.features(), self.read_length());
            builder = builder.set_cigar(cigar);
        }

        if let Some(mate_reference_sequence_name) = get_reference_sequence_name(
            reference_sequences,
            self.next_fragment_reference_sequence_id(),
        )? {
            builder = builder.set_mate_reference_sequence_name(mate_reference_sequence_name);
        }

        if let Some(mate_alignment_start) = self.mate_alignment_start() {
            let mate_alignment_start =
                sam::record::Position::try_from(usize::from(mate_alignment_start) as i32).unwrap();
            builder = builder.set_mate_position(mate_alignment_start);
        }

        builder = builder.set_template_length(self.template_size());

        if !self.bases().is_empty() {
            builder = builder.set_sequence(self.bases().clone());
        }

        if !self.quality_scores().is_empty() {
            builder = builder.set_quality_scores(self.quality_scores().clone());
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
