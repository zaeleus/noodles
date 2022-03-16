use std::io;

use noodles_bam as bam;
use noodles_sam::{self as sam, AlignmentRecord};

use super::{resolve::resolve_features, Features, Record};

impl Record {
    /// Converts an alignment record to a CRAM record.
    pub fn try_from_alignment_record<R>(header: &sam::Header, record: &R) -> io::Result<Self>
    where
        R: AlignmentRecord,
    {
        let mut builder = Self::builder();

        builder = builder.set_bam_flags(record.flags());

        // CRAM flags

        if let Some(reference_sequence) = record
            .reference_sequence(header.reference_sequences())
            .transpose()?
        {
            let reference_sequence_id =
                get_reference_sequence_id(header.reference_sequences(), reference_sequence)?;
            builder = builder.set_reference_sequence_id(reference_sequence_id);
        }

        builder = builder.set_read_length(record.sequence().len());

        if let Some(alignment_start) = record.alignment_start() {
            builder = builder.set_alignment_start(alignment_start);
        }

        if let Some(read_group_id) = get_read_group_id(header.read_groups(), record.data())? {
            builder = builder.set_read_group_id(read_group_id);
        }

        if let Some(read_name) = record.read_name() {
            builder = builder.set_read_name(read_name.clone());
        }

        // next mate bit flags

        if let Some(reference_sequence) = record
            .mate_reference_sequence(header.reference_sequences())
            .transpose()?
        {
            let reference_sequence_id =
                get_reference_sequence_id(header.reference_sequences(), reference_sequence)?;
            builder = builder.set_next_fragment_reference_sequence_id(reference_sequence_id);
        }

        if let Some(mate_alignment_start) = record.mate_alignment_start() {
            builder = builder.set_next_mate_alignment_start(mate_alignment_start);
        }

        builder = builder.set_template_size(record.template_length());

        // distance to next fragment

        if !record.data().is_empty() {
            use sam::record::data::field::Tag;
            let mut data = record.data().clone();
            data.remove(Tag::ReadGroup);
            builder = builder.set_tags(data);
        }

        builder = builder.set_bases(record.sequence().clone());

        let features = Features::from_cigar(record.cigar(), record.sequence());
        builder = builder.set_features(features);

        if let Some(mapping_quality) = record.mapping_quality() {
            builder = builder.set_mapping_quality(mapping_quality);
        }

        builder = builder.set_quality_scores(record.quality_scores().clone());

        Ok(builder.build())
    }

    /// Converts this CRAM record to a SAM record.
    ///
    /// This assumes this record is fully resolved.
    pub fn try_into_sam_record(&self, header: &sam::Header) -> io::Result<sam::Record> {
        let mut builder = sam::Record::builder();

        if let Some(read_name) = self.read_name() {
            builder = builder.set_read_name(read_name.clone());
        }

        builder = builder.set_flags(self.bam_flags());

        if let Some(reference_sequence_name) =
            get_reference_sequence_name(header.reference_sequences(), self.reference_sequence_id())?
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
            let cigar = resolve_features(self.features(), self.read_length());
            builder = builder.set_cigar(cigar);
        }

        if let Some(mate_reference_sequence_name) = get_reference_sequence_name(
            header.reference_sequences(),
            self.next_fragment_reference_sequence_id(),
        )? {
            builder = builder.set_mate_reference_sequence_name(mate_reference_sequence_name);
        }

        if let Some(mate_alignment_start) = self.mate_alignment_start() {
            builder = builder.set_mate_position(mate_alignment_start);
        }

        builder = builder.set_template_length(self.template_size());

        if !self.bases().is_empty() {
            builder = builder.set_sequence(self.bases().clone());
        }

        if !self.quality_scores().is_empty() {
            builder = builder.set_quality_scores(self.quality_scores().clone());
        }

        let data = build_data(header.read_groups(), self.tags(), self.read_group_id())?;
        builder = builder.set_data(data);

        Ok(builder.build())
    }
}

fn get_reference_sequence_id(
    reference_sequences: &sam::header::ReferenceSequences,
    reference_sequence: &sam::header::ReferenceSequence,
) -> io::Result<bam::record::ReferenceSequenceId> {
    reference_sequences
        .get_index_of(reference_sequence.name().as_str())
        .map(bam::record::ReferenceSequenceId::from)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid reference sequence name",
            )
        })
}

fn get_read_group_id(
    read_groups: &sam::header::ReadGroups,
    data: &sam::record::Data,
) -> io::Result<Option<usize>> {
    use sam::record::data::field::Tag;

    let rg = match data.get(Tag::ReadGroup) {
        Some(field) => field,
        None => return Ok(None),
    };

    let read_group_name = rg.value().as_str().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid read group field value",
        )
    })?;

    read_groups
        .get_index_of(read_group_name)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid read group name"))
        .map(Some)
}

fn get_reference_sequence_name(
    reference_sequences: &sam::header::ReferenceSequences,
    reference_sequence_id: Option<bam::record::ReferenceSequenceId>,
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

fn build_data(
    read_groups: &sam::header::ReadGroups,
    tags: &sam::record::Data,
    read_group_id: Option<usize>,
) -> io::Result<sam::record::Data> {
    use sam::record::data::{
        field::{Tag as SamTag, Value},
        Field,
    };

    let mut data = tags.clone();

    if let Some(id) = read_group_id {
        let name = read_groups
            .get_index(id)
            .map(|(_, rg)| rg.id())
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid read group ID"))?;

        let rg = Field::new(SamTag::ReadGroup, Value::String(name.into()));
        data.insert(rg);
    }

    Ok(data)
}
