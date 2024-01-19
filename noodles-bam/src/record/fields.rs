//! BAM record fields.

mod cigar;
pub mod data;
mod flags;
mod mapping_quality;
mod name;
mod position;
mod quality_scores;
mod reference_sequence_id;
mod sequence;
mod template_length;

use std::mem;

pub use self::{
    cigar::Cigar, data::Data, flags::Flags, mapping_quality::MappingQuality, name::Name,
    position::Position, quality_scores::QualityScores, reference_sequence_id::ReferenceSequenceId,
    sequence::Sequence, template_length::TemplateLength,
};
use super::{bounds, Bounds};

pub(super) struct Fields<'a> {
    buf: &'a [u8],
    bounds: &'a Bounds,
}

impl<'a> Fields<'a> {
    pub(super) fn new(buf: &'a [u8], bounds: &'a Bounds) -> Self {
        Self { buf, bounds }
    }

    pub(super) fn reference_sequence_id(&self) -> Option<ReferenceSequenceId> {
        let src = &self.buf[bounds::REFERENCE_SEQUENCE_ID_RANGE];
        // SAFETY: `src` is 4 bytes.
        get_reference_sequence_id(src.try_into().unwrap())
    }

    pub(super) fn alignment_start(&self) -> Option<Position> {
        let src = &self.buf[bounds::ALIGNMENT_START_RANGE];
        // SAFETY: `src` is 4 bytes.
        get_position(src.try_into().unwrap())
    }

    pub(super) fn mapping_quality(&self) -> Option<MappingQuality> {
        const MISSING: u8 = 255;

        match self.buf[bounds::MAPPING_QUALITY_INDEX] {
            MISSING => None,
            n => Some(MappingQuality::new(n)),
        }
    }

    pub(super) fn flags(&self) -> Flags {
        let src = &self.buf[bounds::FLAGS_RANGE];
        // SAFETY: `src` is 2 bytes.
        let n = u16::from_le_bytes(src.try_into().unwrap());
        Flags::new(n)
    }

    pub(super) fn mate_reference_sequence_id(&self) -> Option<ReferenceSequenceId> {
        let src = &self.buf[bounds::MATE_REFERENCE_SEQUENCE_ID_RANGE];
        // SAFETY: `src` is 4 bytes.
        get_reference_sequence_id(src.try_into().unwrap())
    }

    pub(super) fn mate_alignment_start(&self) -> Option<Position> {
        let src = &self.buf[bounds::MATE_ALIGNMENT_START_RANGE];
        get_position(src.try_into().unwrap())
    }

    pub(super) fn template_length(&self) -> TemplateLength {
        let src = &self.buf[bounds::TEMPLATE_LENGTH_RANGE];
        // SAFETY: `src` is 4 bytes.
        let n = i32::from_le_bytes(src.try_into().unwrap());
        TemplateLength::new(n)
    }

    pub(super) fn name(&self) -> Option<Name<'a>> {
        const MISSING: &[u8] = &[b'*', 0x00];

        match &self.buf[self.bounds.name_range()] {
            MISSING => None,
            buf => Some(Name::new(buf)),
        }
    }

    pub(super) fn cigar(&self) -> Cigar<'a> {
        use bytes::Buf;

        use self::data::get_raw_cigar;

        const SKIP: u8 = 3;
        const SOFT_CLIP: u8 = 4;

        fn decode_op(n: u32) -> (u8, usize) {
            ((n & 0x0f) as u8, usize::try_from(n >> 4).unwrap())
        }

        let mut src = &self.buf[self.bounds.cigar_range()];

        if src.len() == 2 * mem::size_of::<u32>() {
            let k = self.sequence().len();

            // SAFETY: `src` is 8 bytes.
            let op_1 = decode_op(src.get_u32_le());
            let op_2 = decode_op(src.get_u32_le());

            if op_1 == (SOFT_CLIP, k) && matches!(op_2, (SKIP, _)) {
                let mut data_src = &self.buf[self.bounds.data_range()];

                if let Ok(Some(buf)) = get_raw_cigar(&mut data_src) {
                    return Cigar::new(buf);
                }
            }
        }

        Cigar::new(src)
    }

    pub(super) fn sequence(&self) -> Sequence<'a> {
        let src = &self.buf[self.bounds.sequence_range()];
        let quality_scores_range = self.bounds.quality_scores_range();
        let base_count = quality_scores_range.end - quality_scores_range.start;
        Sequence::new(src, base_count)
    }

    pub(super) fn quality_scores(&self) -> QualityScores<'a> {
        let src = &self.buf[self.bounds.quality_scores_range()];
        QualityScores::new(src)
    }

    pub(super) fn data(&self) -> Data<'a> {
        let src = &self.buf[self.bounds.data_range()];
        Data::new(src)
    }
}

fn get_reference_sequence_id(src: [u8; 4]) -> Option<ReferenceSequenceId> {
    const UNMAPPED: i32 = -1;

    match i32::from_le_bytes(src) {
        UNMAPPED => None,
        n => Some(ReferenceSequenceId::new(n)),
    }
}

fn get_position(src: [u8; 4]) -> Option<Position> {
    const MISSING: i32 = -1;

    match i32::from_le_bytes(src) {
        MISSING => None,
        n => Some(Position::new(n)),
    }
}