use crate::{Cigar, Data, Flags, MappingQuality};

use super::{Record, NULL_FIELD};

#[derive(Debug, Default)]
pub struct Builder {
    name: Option<String>,
    flags: Flags,
    reference_sequence_name: Option<String>,
    position: u32,
    mapping_quality: MappingQuality,
    cigar: Option<Cigar>,
    mate_reference_sequence_name: Option<String>,
    mate_position: u32,
    template_len: i32,
    sequence: Option<String>,
    quality_scores: Option<String>,
    data: Option<Data>,
}

impl Builder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_name(mut self, name: &str) -> Self {
        self.name = Some(name.into());
        self
    }

    pub fn set_flags(mut self, flags: Flags) -> Self {
        self.flags = flags;
        self
    }

    pub fn set_reference_sequence_name(mut self, reference_sequence_name: &str) -> Self {
        self.reference_sequence_name = Some(reference_sequence_name.into());
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
        self.cigar = Some(cigar);
        self
    }

    pub fn set_mate_reference_sequence_name(mut self, mate_reference_sequence_name: &str) -> Self {
        self.mate_reference_sequence_name = Some(mate_reference_sequence_name.into());
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

    pub fn set_sequence(mut self, sequence: &str) -> Self {
        self.sequence = Some(sequence.into());
        self
    }

    pub fn set_quality_scores(mut self, quality_scores: &str) -> Self {
        self.quality_scores = Some(quality_scores.into());
        self
    }

    pub fn set_data(mut self, data: Data) -> Self {
        self.data = Some(data);
        self
    }

    pub fn build(self) -> Record {
        let null_field = || NULL_FIELD.into();

        Record {
            qname: self.name.unwrap_or_else(null_field),
            flag: self.flags,
            rname: self.reference_sequence_name.unwrap_or_else(null_field),
            pos: self.position,
            mapq: self.mapping_quality,
            cigar: self.cigar.unwrap_or_default(),
            rnext: self.mate_reference_sequence_name.unwrap_or_else(null_field),
            pnext: self.mate_position,
            tlen: self.template_len,
            seq: self.sequence.unwrap_or_else(null_field),
            qual: self.quality_scores.unwrap_or_else(null_field),
            data: self.data.unwrap_or_default(),
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

        assert_eq!(record.name(), "*");
        assert_eq!(u16::from(record.flags()), 0);
        assert_eq!(record.reference_sequence_name(), "*");
        assert_eq!(record.position(), 0);
        assert_eq!(u8::from(record.mapping_quality()), 255);
        assert!(record.cigar().ops().is_empty());
        assert_eq!(record.mate_reference_sequence_name(), "*");
        assert_eq!(record.mate_position(), 0);
        assert_eq!(record.template_len(), 0);
        assert_eq!(record.sequence(), "*");
        assert_eq!(record.quality_scores(), "*");
        assert!(record.data().fields().is_empty());
    }

    #[test]
    fn test_build() {
        let cigar = Cigar::new(vec![cigar::Op::new(cigar::op::Kind::Match, 4)]);

        let data = Data::new(vec![data::Field::new(
            String::from("NH"),
            data::Value::Int32(1),
        )]);

        let record = Builder::new()
            .set_name("r0")
            .set_flags(Flags::from(65))
            .set_reference_sequence_name("sq0")
            .set_position(13)
            .set_mapping_quality(MappingQuality::from(37))
            .set_cigar(cigar)
            .set_mate_reference_sequence_name("sq1")
            .set_mate_position(17)
            .set_template_len(4)
            .set_sequence("ATCGATC")
            .set_quality_scores("NOODLES")
            .set_data(data)
            .build();

        assert_eq!(record.name(), "r0");
        assert_eq!(u16::from(record.flags()), 65);
        assert_eq!(record.reference_sequence_name(), "sq0");
        assert_eq!(record.position(), 13);
        assert_eq!(u8::from(record.mapping_quality()), 37);
        assert_eq!(record.cigar().ops().len(), 1);
        assert_eq!(record.mate_reference_sequence_name(), "sq1");
        assert_eq!(record.mate_position(), 17);
        assert_eq!(record.template_len(), 4);
        assert_eq!(record.sequence(), "ATCGATC");
        assert_eq!(record.quality_scores(), "NOODLES");
        assert_eq!(record.data().fields().len(), 1);
    }
}
