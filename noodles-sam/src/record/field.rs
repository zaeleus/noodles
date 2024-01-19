//! SAM record field.

mod cigar;
mod data;
mod flags;
mod mapping_quality;
mod name;
mod position;
mod quality_scores;
mod reference_sequence_id;
mod reference_sequence_name;
mod sequence;
mod template_length;

pub use self::{
    cigar::Cigar, data::Data, flags::Flags, mapping_quality::MappingQuality, name::Name,
    position::Position, quality_scores::QualityScores, reference_sequence_id::ReferenceSequenceId,
    reference_sequence_name::ReferenceSequenceName, sequence::Sequence,
    template_length::TemplateLength,
};
use super::Bounds;
use crate::Header;

const MISSING: &[u8] = b"*";

pub(super) struct Fields<'a> {
    buf: &'a [u8],
    bounds: &'a Bounds,
}

impl<'a> Fields<'a> {
    pub(super) fn new(buf: &'a [u8], bounds: &'a Bounds) -> Self {
        Self { buf, bounds }
    }

    pub fn name(&self) -> Option<Name<'a>> {
        match &self.buf[self.bounds.name_range()] {
            MISSING => None,
            buf => Some(Name::new(buf)),
        }
    }

    pub fn flags(&self) -> Flags<'a> {
        let src = &self.buf[self.bounds.flags_range()];
        Flags::new(src)
    }

    pub fn reference_sequence_id<'h: 'a>(
        &self,
        header: &'h Header,
    ) -> Option<ReferenceSequenceId<'h, 'a>> {
        self.reference_sequence_name()
            .map(|reference_sequence_name| {
                ReferenceSequenceId::new(header, reference_sequence_name)
            })
    }

    pub fn reference_sequence_name(&self) -> Option<ReferenceSequenceName<'a>> {
        match &self.buf[self.bounds.reference_sequence_name_range()] {
            MISSING => None,
            buf => Some(ReferenceSequenceName::new(buf)),
        }
    }

    pub fn alignment_start(&self) -> Option<Position<'a>> {
        const MISSING: &[u8] = b"0";

        match &self.buf[self.bounds.alignment_start_range()] {
            MISSING => None,
            buf => Some(Position::new(buf)),
        }
    }

    pub fn mapping_quality(&self) -> Option<MappingQuality<'a>> {
        const MISSING: &[u8] = b"255";

        match &self.buf[self.bounds.mapping_quality_range()] {
            MISSING => None,
            buf => Some(MappingQuality::new(buf)),
        }
    }

    pub fn cigar(&self) -> Cigar<'a> {
        match &self.buf[self.bounds.cigar_range()] {
            MISSING => Cigar::new(b""),
            buf => Cigar::new(buf),
        }
    }

    pub fn mate_reference_sequence_id<'h: 'a>(
        &self,
        header: &'h Header,
    ) -> Option<ReferenceSequenceId<'h, 'a>> {
        self.mate_reference_sequence_name()
            .map(|mate_reference_sequence_name| {
                ReferenceSequenceId::new(header, mate_reference_sequence_name)
            })
    }

    pub fn mate_reference_sequence_name(&self) -> Option<ReferenceSequenceName<'a>> {
        const EQ: &[u8] = b"=";

        match &self.buf[self.bounds.mate_reference_sequence_name_range()] {
            MISSING => None,
            EQ => self.reference_sequence_name(),
            buf => Some(ReferenceSequenceName::new(buf)),
        }
    }

    pub fn mate_alignment_start(&self) -> Option<Position<'a>> {
        const MISSING: &[u8] = b"0";

        match &self.buf[self.bounds.mate_alignment_start_range()] {
            MISSING => None,
            buf => Some(Position::new(buf)),
        }
    }

    pub fn template_length(&self) -> TemplateLength<'a> {
        let buf = &self.buf[self.bounds.template_length_range()];
        TemplateLength::new(buf)
    }

    pub fn sequence(&self) -> Sequence<'a> {
        let buf = match &self.buf[self.bounds.sequence_range()] {
            MISSING => b"",
            buf => buf,
        };

        Sequence::new(buf)
    }

    pub fn quality_scores(&self) -> QualityScores<'a> {
        let buf = match &self.buf[self.bounds.quality_scores_range()] {
            MISSING => b"",
            buf => buf,
        };

        QualityScores::new(buf)
    }

    pub fn data(&self) -> Data<'a> {
        let buf = &self.buf[self.bounds.data_range()];
        Data::new(buf)
    }
}
