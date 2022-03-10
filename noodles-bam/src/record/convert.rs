use std::io;

use noodles_sam::{self as sam, AlignmentRecord};

use super::{Record, ReferenceSequenceId};

impl Record {
    /// Converts a SAM record to a BAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam as bam;
    /// use noodles_sam::{self as sam, header::ReferenceSequences};
    ///
    /// let reference_sequences = ReferenceSequences::default();
    /// let sam_record = sam::Record::default();
    ///
    /// let record = bam::Record::try_from_sam_record(&reference_sequences, &sam_record)?;
    /// assert_eq!(record, bam::Record::default());
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn try_from_sam_record(
        reference_sequences: &sam::header::ReferenceSequences,
        sam_record: &sam::Record,
    ) -> io::Result<Self> {
        use crate::{reader::record::read_record, writer::sam_record::write_sam_record};

        let mut buf = Vec::new();
        write_sam_record(&mut buf, reference_sequences, sam_record)?;

        let mut reader = &buf[..];
        let mut record = Self::default();
        read_record(&mut reader, &mut Vec::new(), &mut record)?;

        Ok(record)
    }

    /// Converts this record to a SAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam as bam;
    /// use noodles_sam as sam;
    ///
    /// let reference_sequences = sam::header::ReferenceSequences::default();
    ///
    /// let record = bam::Record::default();
    /// let sam_record = record.try_into_sam_record(&reference_sequences)?;
    ///
    /// assert_eq!(sam_record, sam::Record::default());
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn try_into_sam_record(
        &self,
        reference_sequences: &sam::header::ReferenceSequences,
    ) -> io::Result<sam::Record> {
        let mut builder = sam::Record::builder();

        if let Some(read_name) = self.read_name() {
            builder = builder.set_read_name(read_name.clone());
        }

        builder = builder.set_flags(self.flags());

        if let Some(reference_sequence) =
            self.reference_sequence(reference_sequences).transpose()?
        {
            builder = builder.set_reference_sequence_name(reference_sequence.name().clone());
        }

        if let Some(position) = self.position() {
            let position = sam::record::Position::try_from(usize::from(position) as i32).unwrap();
            builder = builder.set_position(position);
        }

        if let Some(mapping_quality) = self.mapping_quality() {
            builder = builder.set_mapping_quality(mapping_quality);
        }

        if !self.cigar().is_empty() {
            builder = builder.set_cigar(self.cigar().clone());
        }

        if let Some(mate_reference_sequence_name) =
            get_reference_sequence_name(reference_sequences, self.mate_reference_sequence_id())?
        {
            builder = builder.set_mate_reference_sequence_name(mate_reference_sequence_name);
        }

        if let Some(mate_position) = self.mate_position() {
            builder = builder.set_mate_position(mate_position);
        }

        builder = builder.set_template_length(self.template_length());

        if !self.sequence().is_empty() {
            builder = builder.set_sequence(self.sequence().clone());
        }

        if !self.quality_scores().is_empty() {
            builder = builder.set_quality_scores(self.quality_scores().clone());
        }

        let data = self
            .data()
            .try_into()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        builder = builder.set_data(data);

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

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;

    fn build_reference_sequences(
    ) -> Result<sam::header::ReferenceSequences, Box<dyn std::error::Error>> {
        use sam::header::{reference_sequence, ReferenceSequence};

        [
            ("sq0".parse()?, 5),
            ("sq1".parse()?, 8),
            ("sq2".parse()?, 13),
        ]
        .into_iter()
        .map(|(name, len): (reference_sequence::Name, i32)| {
            let sn = name.to_string();
            ReferenceSequence::new(name, len).map(|rs| (sn, rs))
        })
        .collect::<Result<_, _>>()
        .map_err(|e| e.into())
    }

    fn build_record() -> Result<Record, Box<dyn std::error::Error>> {
        use sam::record::{Flags, MappingQuality};

        use crate::record::Data;

        let reference_sequence_id = ReferenceSequenceId::from(1);

        let record = Record::builder()
            .set_reference_sequence_id(reference_sequence_id)
            .set_position(Position::try_from(61062)?)
            .set_mapping_quality(MappingQuality::try_from(12)?)
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .set_mate_reference_sequence_id(reference_sequence_id)
            .set_mate_position(sam::record::Position::try_from(61153)?)
            .set_template_length(166)
            .set_read_name("r0".parse()?)
            .set_cigar("4M".parse()?)
            .set_sequence("ATGC".parse()?)
            .set_quality_scores("@>?A".parse()?)
            .set_data(Data::try_from(vec![
                0x4e, 0x4d, 0x43, 0x00, // NM:i:0
                0x50, 0x47, 0x5a, 0x53, 0x4e, 0x41, 0x50, 0x00, // PG:Z:SNAP
            ])?)
            .build()?;

        Ok(record)
    }

    #[test]
    fn test_try_into_sam_record() -> Result<(), Box<dyn std::error::Error>> {
        use sam::record::data::{
            field::{Tag, Value},
            Field,
        };

        let bam_record = build_record()?;
        let reference_sequences = build_reference_sequences()?;
        let actual = bam_record.try_into_sam_record(&reference_sequences)?;

        let expected = sam::Record::builder()
            .set_read_name("r0".parse()?)
            .set_flags(sam::record::Flags::SEGMENTED | sam::record::Flags::FIRST_SEGMENT)
            .set_reference_sequence_name("sq1".parse()?)
            .set_position(sam::record::Position::try_from(61062)?)
            .set_mapping_quality(sam::record::MappingQuality::try_from(12)?)
            .set_cigar("4M".parse()?)
            .set_mate_reference_sequence_name("sq1".parse()?)
            .set_mate_position(sam::record::Position::try_from(61153)?)
            .set_template_length(166)
            .set_sequence("ATGC".parse()?)
            .set_quality_scores("@>?A".parse()?)
            .set_data(sam::record::Data::try_from(vec![
                Field::new(Tag::EditDistance, Value::Int(0)),
                Field::new(Tag::Program, Value::String(String::from("SNAP"))),
            ])?)
            .build()?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
