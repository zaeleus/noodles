//! Alignment record.

pub mod cigar;
pub mod data;
mod flags;
mod mapping_quality;
mod name;
mod position;
mod quality_scores;
mod reference_sequence_id;
mod sequence;
mod template_length;

use std::io;

use noodles_core as core;

pub use self::{
    cigar::Cigar, data::Data, flags::Flags, mapping_quality::MappingQuality, name::Name,
    position::Position, quality_scores::QualityScores, reference_sequence_id::ReferenceSequenceId,
    sequence::Sequence, template_length::TemplateLength,
};
use crate::{
    header::{
        record::value::{map::ReferenceSequence, Map},
        ReferenceSequences,
    },
    Header,
};

/// An alignment record.
pub trait Record {
    /// Returns the name.
    fn name(&self) -> Option<Box<dyn Name + '_>>;

    /// Returns the flags.
    fn flags(&self) -> Box<dyn Flags + '_>;

    /// Returns the reference sequence ID.
    fn reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        header: &'h Header,
    ) -> Option<Box<dyn ReferenceSequenceId + 'r>>;

    /// Returns the alignment start.
    fn alignment_start(&self) -> Option<Box<dyn Position + '_>>;

    /// Returns the mapping quality.
    fn mapping_quality(&self) -> Option<Box<dyn MappingQuality + '_>>;

    /// Returns the CIGAR operations.
    fn cigar(&self) -> Box<dyn Cigar + '_>;

    /// Returns the mate reference sequence ID.
    fn mate_reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        header: &'h Header,
    ) -> Option<Box<dyn ReferenceSequenceId + 'r>>;

    /// Returns the mate alignment start.
    fn mate_alignment_start(&self) -> Option<Box<dyn Position + '_>>;

    /// Returns the template length.
    fn template_length(&self) -> Box<dyn TemplateLength + '_>;

    /// Returns the sequence.
    fn sequence(&self) -> Box<dyn Sequence + '_>;

    /// Returns the quality scores.
    fn quality_scores(&self) -> Box<dyn QualityScores + '_>;

    /// Returns the data.
    fn data(&self) -> Box<dyn Data + '_>;

    /// Returns the associated reference sequence.
    fn reference_sequence<'h>(
        &self,
        header: &'h Header,
    ) -> Option<io::Result<(&'h [u8], &'h Map<ReferenceSequence>)>> {
        get_reference_sequence(
            header.reference_sequences(),
            self.reference_sequence_id(header),
        )
    }

    /// Returns the associated mate reference sequence.
    fn mate_reference_sequence<'h>(
        &self,
        header: &'h Header,
    ) -> Option<io::Result<(&'h [u8], &'h Map<ReferenceSequence>)>> {
        get_reference_sequence(
            header.reference_sequences(),
            self.mate_reference_sequence_id(header),
        )
    }

    /// Returns the alignment span.
    fn alignment_span(&self) -> io::Result<usize> {
        self.cigar().alignment_span()
    }

    /// Calculates the end position.
    fn alignment_end(&self) -> Option<io::Result<core::Position>> {
        let alignment_start = self.alignment_start()?;

        let start = match core::Position::try_from(alignment_start.as_ref()) {
            Ok(position) => position,
            Err(e) => return Some(Err(e)),
        };

        let span = match self.alignment_span() {
            Ok(span) => span,
            Err(e) => return Some(Err(e)),
        };

        let end = usize::from(start) + span - 1;
        core::Position::new(end).map(Ok)
    }
}

impl Record for Box<dyn Record> {
    fn name(&self) -> Option<Box<dyn Name + '_>> {
        (**self).name()
    }

    fn flags(&self) -> Box<dyn Flags + '_> {
        (**self).flags()
    }

    fn reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        header: &'h Header,
    ) -> Option<Box<dyn ReferenceSequenceId + 'r>> {
        (**self).reference_sequence_id(header)
    }

    fn alignment_start(&self) -> Option<Box<dyn Position + '_>> {
        (**self).alignment_start()
    }

    fn mapping_quality(&self) -> Option<Box<dyn MappingQuality + '_>> {
        (**self).mapping_quality()
    }

    fn cigar(&self) -> Box<dyn Cigar + '_> {
        (**self).cigar()
    }

    fn mate_reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        header: &'h Header,
    ) -> Option<Box<dyn ReferenceSequenceId + 'r>> {
        (**self).mate_reference_sequence_id(header)
    }

    fn mate_alignment_start(&self) -> Option<Box<dyn Position + '_>> {
        (**self).mate_alignment_start()
    }

    fn template_length(&self) -> Box<dyn TemplateLength + '_> {
        (**self).template_length()
    }

    fn sequence(&self) -> Box<dyn Sequence + '_> {
        (**self).sequence()
    }

    fn quality_scores(&self) -> Box<dyn QualityScores + '_> {
        (**self).quality_scores()
    }

    fn data(&self) -> Box<dyn Data + '_> {
        (**self).data()
    }
}

fn get_reference_sequence<'r, 'h: 'r>(
    reference_sequences: &'h ReferenceSequences,
    reference_sequence_id: Option<Box<dyn ReferenceSequenceId + 'r>>,
) -> Option<io::Result<(&'h [u8], &'h Map<ReferenceSequence>)>> {
    let reference_sequence_id = reference_sequence_id?;

    let id = match reference_sequence_id.try_to_usize() {
        Ok(id) => id,
        Err(e) => return Some(Err(e)),
    };

    let result = reference_sequences
        .get_index(id)
        .map(|(name, reference_sequence)| (name.as_ref(), reference_sequence))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid reference sequence ID"));

    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alignment_end() -> Result<(), Box<dyn std::error::Error>> {
        use crate::alignment::{
            record::cigar::{op::Kind, Op},
            RecordBuf,
        };

        let record = RecordBuf::builder()
            .set_alignment_start(core::Position::try_from(8)?)
            .set_cigar([Op::new(Kind::Match, 5)].into_iter().collect())
            .build();

        let actual = Record::alignment_end(&record).transpose()?;
        let expected = core::Position::new(12);
        assert_eq!(actual, expected);

        Ok(())
    }
}
