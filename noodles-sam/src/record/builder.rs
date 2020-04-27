use crate::{Cigar, Data, Flags, MappingQuality};

use super::{
    MateReferenceSequenceName, QualityScores, ReadName, Record, ReferenceSequenceName, Sequence,
};

#[derive(Debug, Default)]
pub struct Builder {
    name: ReadName,
    flags: Flags,
    reference_sequence_name: ReferenceSequenceName,
    position: u32,
    mapping_quality: MappingQuality,
    cigar: Cigar,
    mate_reference_sequence_name: MateReferenceSequenceName,
    mate_position: u32,
    template_len: i32,
    sequence: Sequence,
    quality_scores: QualityScores,
    data: Data,
}

impl Builder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_name(mut self, name: ReadName) -> Self {
        self.name = name;
        self
    }

    pub fn set_flags(mut self, flags: Flags) -> Self {
        self.flags = flags;
        self
    }

    pub fn set_reference_sequence_name(
        mut self,
        reference_sequence_name: ReferenceSequenceName,
    ) -> Self {
        self.reference_sequence_name = reference_sequence_name;
        self
    }

    pub fn set_position(mut self, position: u32) -> Self {
        self.position = position;
        self
    }

    pub fn set_mapping_quality(mut self, mapping_quality: MappingQuality) -> Self {
        self.mapping_quality = mapping_quality;
        self
    }

    pub fn set_cigar(mut self, cigar: Cigar) -> Self {
        self.cigar = cigar;
        self
    }

    pub fn set_mate_reference_sequence_name(
        mut self,
        mate_reference_sequence_name: MateReferenceSequenceName,
    ) -> Self {
        self.mate_reference_sequence_name = mate_reference_sequence_name;
        self
    }

    pub fn set_mate_position(mut self, mate_position: u32) -> Self {
        self.mate_position = mate_position;
        self
    }

    pub fn set_template_len(mut self, template_len: i32) -> Self {
        self.template_len = template_len;
        self
    }

    pub fn set_sequence(mut self, sequence: Sequence) -> Self {
        self.sequence = sequence;
        self
    }

    pub fn set_quality_scores(mut self, quality_scores: QualityScores) -> Self {
        self.quality_scores = quality_scores;
        self
    }

    pub fn set_data(mut self, data: Data) -> Self {
        self.data = data;
        self
    }

    pub fn build(self) -> Record {
        Record {
            qname: self.name,
            flag: self.flags,
            rname: self.reference_sequence_name,
            pos: self.position,
            mapq: self.mapping_quality,
            cigar: self.cigar,
            rnext: self.mate_reference_sequence_name,
            pnext: self.mate_position,
            tlen: self.template_len,
            seq: self.sequence,
            qual: self.quality_scores,
            data: self.data,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{cigar, data};

    use super::*;

    #[test]
    fn test_default() {
        let record = Builder::new().build();

        assert!(record.name().is_none());
        assert!(record.flags().is_empty());
        assert!(record.reference_sequence_name().is_none());
        assert_eq!(record.position(), 0);
        assert_eq!(u8::from(record.mapping_quality()), 255);
        assert!(record.cigar().ops().is_empty());
        assert!(record.mate_reference_sequence_name().is_none());
        assert_eq!(record.mate_position(), 0);
        assert_eq!(record.template_len(), 0);
        assert!(record.sequence().is_empty());
        assert!(record.quality_scores().is_empty());
        assert!(record.data().fields().is_empty());
    }

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        let name: ReadName = "r0".parse()?;
        let reference_sequence_name: ReferenceSequenceName = "sq0".parse()?;
        let cigar = Cigar::new(vec![cigar::Op::new(cigar::op::Kind::Match, 4)]);
        let mate_reference_sequence_name: MateReferenceSequenceName = MateReferenceSequenceName::Eq;
        let sequence: Sequence = "ATCGATC".parse()?;
        let quality_scores: QualityScores = "NOODLES".parse()?;

        let data = Data::new(vec![data::Field::new(
            String::from("NH"),
            data::Value::Int32(1),
        )]);

        let record = Builder::new()
            .set_name(name.clone())
            .set_flags(Flags::from(65))
            .set_reference_sequence_name(reference_sequence_name.clone())
            .set_position(13)
            .set_mapping_quality(MappingQuality::from(37))
            .set_cigar(cigar)
            .set_mate_reference_sequence_name(mate_reference_sequence_name.clone())
            .set_mate_position(17)
            .set_template_len(4)
            .set_sequence(sequence.clone())
            .set_quality_scores(quality_scores.clone())
            .set_data(data)
            .build();

        assert_eq!(record.name(), &name);
        assert_eq!(u16::from(record.flags()), 65);
        assert_eq!(record.reference_sequence_name(), &reference_sequence_name);
        assert_eq!(record.position(), 13);
        assert_eq!(u8::from(record.mapping_quality()), 37);
        assert_eq!(record.cigar().ops().len(), 1);

        assert_eq!(
            record.mate_reference_sequence_name(),
            &mate_reference_sequence_name
        );

        assert_eq!(record.mate_position(), 17);
        assert_eq!(record.template_len(), 4);
        assert_eq!(record.sequence(), &sequence);
        assert_eq!(record.quality_scores(), &quality_scores);
        assert_eq!(record.data().fields().len(), 1);

        Ok(())
    }
}
