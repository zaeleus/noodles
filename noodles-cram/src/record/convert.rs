use std::io;

use noodles_sam::{
    self as sam,
    header::{record::value::map, ReferenceSequences},
};

use super::{Features, Flags, QualityScores, Record, Sequence};

impl Record {
    /// Converts an alignment record to a CRAM record.
    pub fn try_from_alignment_record(
        header: &sam::Header,
        record: &sam::alignment::RecordBuf,
    ) -> io::Result<Self> {
        let mut builder = Self::builder();

        let bam_flags = record.flags();
        builder = builder.set_bam_flags(bam_flags);

        let mut flags = Flags::default();

        if let Some((reference_sequence_name, _)) = record.reference_sequence(header).transpose()? {
            let reference_sequence_id =
                get_reference_sequence_id(header.reference_sequences(), reference_sequence_name)?;
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

        if let Some((reference_sequence_name, _)) =
            record.mate_reference_sequence(header).transpose()?
        {
            let reference_sequence_id =
                get_reference_sequence_id(header.reference_sequences(), reference_sequence_name)?;
            builder = builder.set_next_fragment_reference_sequence_id(reference_sequence_id);
        }

        if let Some(mate_alignment_start) = record.mate_alignment_start() {
            builder = builder.set_next_mate_alignment_start(mate_alignment_start);
        }

        builder = builder.set_template_size(record.template_length());

        // distance to next fragment

        if !record.data().is_empty() {
            use sam::record::data::field::tag;
            let mut data = record.data().clone();
            data.remove(&tag::READ_GROUP);
            builder = builder.set_tags(data);
        }

        let raw_bases: Vec<_> = record
            .sequence()
            .as_ref()
            .iter()
            .copied()
            .map(u8::from)
            .collect();
        let bases = Sequence::from(raw_bases);

        builder = builder.set_bases(bases);

        if !bam_flags.is_unmapped() {
            let features = Features::from_cigar(
                flags,
                record.cigar(),
                record.sequence(),
                record.quality_scores(),
            );

            builder = builder.set_features(features);
        }

        if let Some(mapping_quality) = record.mapping_quality() {
            builder = builder.set_mapping_quality(mapping_quality);
        }

        if !record.quality_scores().is_empty() {
            if bam_flags.is_unmapped() {
                flags.insert(Flags::QUALITY_SCORES_STORED_AS_ARRAY);
            }

            let scores: Vec<_> = record
                .quality_scores()
                .as_ref()
                .iter()
                .copied()
                .map(u8::from)
                .collect();

            let quality_scores = QualityScores::from(scores);

            builder = builder.set_quality_scores(quality_scores);
        }

        Ok(builder.set_flags(flags).build())
    }

    /// Converts this CRAM record to an alignment record.
    pub fn try_into_alignment_record(
        self,
        header: &sam::Header,
    ) -> io::Result<sam::alignment::RecordBuf> {
        let mut builder = sam::alignment::RecordBuf::builder();

        if let Some(read_name) = self.read_name {
            builder = builder.set_read_name(read_name);
        }

        builder = builder.set_flags(self.bam_bit_flags);

        if let Some(reference_sequence_id) = self.reference_sequence_id {
            builder = builder.set_reference_sequence_id(reference_sequence_id);
        }

        if let Some(alignment_start) = self.alignment_start {
            builder = builder.set_alignment_start(alignment_start);
        }

        if let Some(mapping_quality) = self.mapping_quality {
            builder = builder.set_mapping_quality(mapping_quality);
        }

        if !self.bam_bit_flags.is_unmapped() {
            let cigar = self
                .features
                .try_into_cigar(self.read_length)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

            builder = builder.set_cigar(cigar);
        }

        if let Some(mate_reference_sequence_id) = self.next_fragment_reference_sequence_id {
            builder = builder.set_mate_reference_sequence_id(mate_reference_sequence_id);
        }

        if let Some(mate_alignment_start) = self.next_mate_alignment_start {
            builder = builder.set_mate_alignment_start(mate_alignment_start);
        }

        builder = builder.set_template_length(self.template_size);

        if !self.bases.is_empty() {
            let raw_bases = Vec::<_>::from(self.bases);
            let bases = sam::record::Sequence::try_from(raw_bases)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
            builder = builder.set_sequence(bases);
        }

        builder = builder.set_quality_scores(self.quality_scores);

        let mut data = self.tags;
        maybe_insert_read_group(&mut data, header.read_groups(), self.read_group_id)?;
        builder = builder.set_data(data);

        Ok(builder.build())
    }
}

fn get_reference_sequence_id(
    reference_sequences: &ReferenceSequences,
    reference_sequence_name: &map::reference_sequence::Name,
) -> io::Result<usize> {
    reference_sequences
        .get_index_of(reference_sequence_name)
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
    use sam::record::data::field::tag;

    let Some(rg_value) = data.get(&tag::READ_GROUP) else {
        return Ok(None);
    };

    let read_group_name = rg_value.as_str().ok_or_else(|| {
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

fn maybe_insert_read_group(
    data: &mut sam::record::Data,
    read_groups: &sam::header::ReadGroups,
    read_group_id: Option<usize>,
) -> io::Result<()> {
    use sam::record::data::field::{tag, Value};

    if let Some(id) = read_group_id {
        let name = read_groups
            .get_index(id)
            .map(|(name, _)| name)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid read group ID"))?;

        data.insert(tag::READ_GROUP, Value::String(name.into()));
    }

    Ok(())
}
