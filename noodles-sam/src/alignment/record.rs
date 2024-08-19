//! Alignment record.

pub mod cigar;
pub mod data;
mod flags;
pub mod mapping_quality;
mod quality_scores;
mod sequence;

use std::io;

use bstr::BStr;
use noodles_core::Position;

pub use self::{
    cigar::Cigar, data::Data, flags::Flags, mapping_quality::MappingQuality,
    quality_scores::QualityScores, sequence::Sequence,
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
    fn name(&self) -> Option<&BStr>;

    /// Returns the flags.
    fn flags(&self) -> io::Result<Flags>;

    /// Returns the reference sequence ID.
    fn reference_sequence_id<'r, 'h: 'r>(&'r self, header: &'h Header)
        -> Option<io::Result<usize>>;

    /// Returns the alignment start.
    ///
    /// This position is 1-based, inclusive.
    fn alignment_start(&self) -> Option<io::Result<Position>>;

    /// Returns the mapping quality.
    fn mapping_quality(&self) -> Option<io::Result<MappingQuality>>;

    /// Returns the CIGAR operations.
    fn cigar(&self) -> Box<dyn Cigar + '_>;

    /// Returns the mate reference sequence ID.
    fn mate_reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        header: &'h Header,
    ) -> Option<io::Result<usize>>;

    /// Returns the mate alignment start.
    ///
    /// This position is 1-based, inclusive.
    fn mate_alignment_start(&self) -> Option<io::Result<Position>>;

    /// Returns the template length.
    fn template_length(&self) -> io::Result<i32>;

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
    ) -> Option<io::Result<(&'h BStr, &'h Map<ReferenceSequence>)>> {
        let reference_sequence_id = match self.reference_sequence_id(header).transpose() {
            Ok(reference_sequence_id) => reference_sequence_id,
            Err(e) => return Some(Err(e)),
        };

        get_reference_sequence(header.reference_sequences(), reference_sequence_id)
    }

    /// Returns the associated mate reference sequence.
    fn mate_reference_sequence<'h>(
        &self,
        header: &'h Header,
    ) -> Option<io::Result<(&'h BStr, &'h Map<ReferenceSequence>)>> {
        let mate_reference_sequence_id = match self.mate_reference_sequence_id(header).transpose() {
            Ok(id) => id,
            Err(e) => return Some(Err(e)),
        };

        get_reference_sequence(header.reference_sequences(), mate_reference_sequence_id)
    }

    /// Returns the alignment span.
    fn alignment_span(&self) -> Option<io::Result<usize>> {
        match self.cigar().alignment_span() {
            Ok(0) => None,
            Ok(span) => Some(Ok(span)),
            Err(e) => Some(Err(e)),
        }
    }

    /// Calculates the end position.
    ///
    /// This position is 1-based, inclusive.
    fn alignment_end(&self) -> Option<io::Result<Position>> {
        let start = match self.alignment_start().transpose() {
            Ok(position) => position?,
            Err(e) => return Some(Err(e)),
        };

        match self.alignment_span() {
            Some(Ok(span)) => {
                let end = usize::from(start) + span - 1;
                Position::new(end).map(Ok)
            }
            Some(Err(e)) => Some(Err(e)),
            None => Some(Ok(start)),
        }
    }
}

impl Record for Box<dyn Record> {
    fn name(&self) -> Option<&BStr> {
        (**self).name()
    }

    fn flags(&self) -> io::Result<Flags> {
        (**self).flags()
    }

    fn reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        header: &'h Header,
    ) -> Option<io::Result<usize>> {
        (**self).reference_sequence_id(header)
    }

    fn alignment_start(&self) -> Option<io::Result<Position>> {
        (**self).alignment_start()
    }

    fn mapping_quality(&self) -> Option<io::Result<MappingQuality>> {
        (**self).mapping_quality()
    }

    fn cigar(&self) -> Box<dyn Cigar + '_> {
        (**self).cigar()
    }

    fn mate_reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        header: &'h Header,
    ) -> Option<io::Result<usize>> {
        (**self).mate_reference_sequence_id(header)
    }

    fn mate_alignment_start(&self) -> Option<io::Result<Position>> {
        (**self).mate_alignment_start()
    }

    fn template_length(&self) -> io::Result<i32> {
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

fn get_reference_sequence(
    reference_sequences: &ReferenceSequences,
    reference_sequence_id: Option<usize>,
) -> Option<io::Result<(&BStr, &Map<ReferenceSequence>)>> {
    let id = reference_sequence_id?;

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
            .set_alignment_start(Position::try_from(8)?)
            .set_cigar([Op::new(Kind::Match, 5)].into_iter().collect())
            .build();

        let actual = Record::alignment_end(&record).transpose()?;
        let expected = Position::new(12);
        assert_eq!(actual, expected);

        Ok(())
    }
}
