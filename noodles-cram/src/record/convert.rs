use std::io;

use noodles_sam as sam;

use super::{features::Cigar, Record};

impl Record<'_> {
    /// Converts this CRAM record to an alignment record.
    pub fn try_into_alignment_record(
        self,
        _header: &sam::Header,
    ) -> io::Result<sam::alignment::RecordBuf> {
        let data = sam::alignment::Record::data(&self)
            .iter()
            .map(|result| result.and_then(|(tag, value)| value.try_into().map(|v| (tag, v))))
            .collect::<io::Result<_>>()?;

        let mut builder = sam::alignment::RecordBuf::builder().set_data(data);

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

        Ok(builder.build())
    }
}
