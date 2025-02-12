use std::io;

use noodles_sam as sam;

use super::{features::Cigar, Record};

impl Record {
    /// Converts this CRAM record to an alignment record.
    pub fn try_into_alignment_record(
        self,
        header: &sam::Header,
    ) -> io::Result<sam::alignment::RecordBuf> {
        let mut builder = sam::alignment::RecordBuf::builder();

        if let Some(read_name) = self.name {
            builder = builder.set_name(read_name);
        }

        builder = builder.set_flags(self.bam_flags);

        if let Some(reference_sequence_id) = self.reference_sequence_id {
            builder = builder.set_reference_sequence_id(reference_sequence_id);
        }

        if let Some(alignment_start) = self.alignment_start {
            builder = builder.set_alignment_start(alignment_start);
        }

        if let Some(mapping_quality) = self.mapping_quality {
            builder = builder.set_mapping_quality(mapping_quality);
        }

        if !self.bam_flags.is_unmapped() {
            let cigar = Cigar::new(&self.features, self.read_length).collect();
            builder = builder.set_cigar(cigar);
        }

        if let Some(mate_reference_sequence_id) = self.mate_reference_sequence_id {
            builder = builder.set_mate_reference_sequence_id(mate_reference_sequence_id);
        }

        if let Some(mate_alignment_start) = self.mate_alignment_start {
            builder = builder.set_mate_alignment_start(mate_alignment_start);
        }

        builder = builder
            .set_template_length(self.template_length)
            .set_sequence(self.sequence)
            .set_quality_scores(self.quality_scores);

        let mut data = self.tags;
        maybe_insert_read_group(&mut data, header.read_groups(), self.read_group_id)?;
        builder = builder.set_data(data);

        Ok(builder.build())
    }
}

fn maybe_insert_read_group(
    data: &mut sam::alignment::record_buf::Data,
    read_groups: &sam::header::ReadGroups,
    read_group_id: Option<usize>,
) -> io::Result<()> {
    use sam::alignment::{record::data::field::Tag, record_buf::data::field::Value};

    if let Some(id) = read_group_id {
        let name = read_groups
            .get_index(id)
            .map(|(name, _)| name)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid read group ID"))?;

        data.insert(Tag::READ_GROUP, Value::String(name.clone()));
    }

    Ok(())
}
