use std::{cmp, io, num::NonZero};

use noodles_core::Position;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Context {
    reference_sequence_id: usize,
    alignment_start: Position,
    alignment_end: Position,
}

impl Context {
    fn new(
        reference_sequence_id: usize,
        alignment_start: Position,
        alignment_end: Position,
    ) -> Self {
        Self {
            reference_sequence_id,
            alignment_start,
            alignment_end,
        }
    }

    pub fn reference_sequence_id(&self) -> usize {
        self.reference_sequence_id
    }

    pub fn alignment_start(&self) -> Position {
        self.alignment_start
    }

    pub fn alignment_span(&self) -> usize {
        usize::from(self.alignment_end) - usize::from(self.alignment_start) + 1
    }

    pub fn alignment_end(&self) -> Position {
        self.alignment_end
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum ReferenceSequenceContext {
    Some(Context),
    #[default]
    None,
    Many,
}

impl ReferenceSequenceContext {
    pub fn some(
        reference_sequence_id: usize,
        alignment_start: Position,
        alignment_end: Position,
    ) -> Self {
        Self::Some(Context::new(
            reference_sequence_id,
            alignment_start,
            alignment_end,
        ))
    }

    pub fn is_many(&self) -> bool {
        matches!(self, Self::Many)
    }

    pub fn update(
        &mut self,
        reference_sequence_id: Option<usize>,
        alignment_start: Option<Position>,
        alignment_end: Option<Position>,
    ) {
        *self = match (*self, reference_sequence_id, alignment_start, alignment_end) {
            (Self::Some(context), Some(record_id), Some(record_start), Some(record_end)) => {
                if record_id == context.reference_sequence_id() {
                    let start = cmp::min(record_start, context.alignment_start());
                    let end = cmp::max(record_end, context.alignment_end());
                    Self::some(record_id, start, end)
                } else {
                    Self::Many
                }
            }
            (Self::Some(..), ..) => Self::Many,
            (Self::None, Some(_), ..) => Self::Many,
            (Self::None, None, ..) => Self::None,
            (Self::Many, ..) => Self::Many,
        }
    }
}

impl TryFrom<(i32, i64, i64)> for ReferenceSequenceContext {
    type Error = io::Error;

    fn try_from(
        (raw_reference_sequence_id, raw_alignment_start, raw_alignment_span): (i32, i64, i64),
    ) -> Result<Self, Self::Error> {
        const UNMAPPED: i32 = -1;
        const MULTIREF: i32 = -2;

        match raw_reference_sequence_id {
            UNMAPPED => Ok(Self::None),
            MULTIREF => Ok(Self::Many),
            _ => {
                let reference_sequence_id = usize::try_from(raw_reference_sequence_id)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                let alignment_start = usize::try_from(raw_alignment_start)
                    .and_then(Position::try_from)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                let alignment_span = usize::try_from(raw_alignment_span)
                    .and_then(NonZero::try_from)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                let alignment_end = alignment_start
                    // SAFETY: `alignment_span > 0`.
                    .checked_add(usize::from(alignment_span) - 1)
                    .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData))?;

                Ok(Self::some(
                    reference_sequence_id,
                    alignment_start,
                    alignment_end,
                ))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alignment_span() -> Result<(), noodles_core::position::TryFromIntError> {
        let context = Context::new(0, Position::try_from(8)?, Position::try_from(13)?);
        assert_eq!(context.alignment_span(), 6);
        Ok(())
    }

    #[test]
    fn test_update() -> Result<(), noodles_core::position::TryFromIntError> {
        let mut context =
            ReferenceSequenceContext::some(0, Position::try_from(8)?, Position::try_from(13)?);
        context.update(Some(0), Position::new(5), Position::new(21));
        assert_eq!(
            context,
            ReferenceSequenceContext::some(0, Position::try_from(5)?, Position::try_from(21)?)
        );

        let mut context =
            ReferenceSequenceContext::some(0, Position::try_from(8)?, Position::try_from(13)?);
        context.update(None, None, None);
        assert_eq!(context, ReferenceSequenceContext::Many);

        let mut context = ReferenceSequenceContext::None;
        context.update(Some(0), Some(Position::MIN), Some(Position::MIN));
        assert_eq!(context, ReferenceSequenceContext::Many);

        let mut context = ReferenceSequenceContext::None;
        context.update(None, None, None);
        assert_eq!(context, ReferenceSequenceContext::None);

        let mut context = ReferenceSequenceContext::Many;
        context.update(Some(0), Some(Position::MIN), Some(Position::MIN));
        assert_eq!(context, ReferenceSequenceContext::Many);

        let mut context = ReferenceSequenceContext::Many;
        context.update(None, None, None);
        assert_eq!(context, ReferenceSequenceContext::Many);

        Ok(())
    }

    #[test]
    fn test_try_from_i32_i64_i64_for_reference_sequence_context()
    -> Result<(), Box<dyn std::error::Error>> {
        assert_eq!(
            ReferenceSequenceContext::try_from((0, 5i64, 8i64))?,
            ReferenceSequenceContext::some(0, Position::try_from(5)?, Position::try_from(12)?)
        );

        assert_eq!(
            ReferenceSequenceContext::try_from((-1, 0i64, 0i64))?,
            ReferenceSequenceContext::None
        );

        assert_eq!(
            ReferenceSequenceContext::try_from((-2, 0i64, 0i64))?,
            ReferenceSequenceContext::Many
        );

        for triplet in [
            (-3, 0i64, 0i64), // invalid reference sequence ID
            (0, -1i64, 1i64), // invalid alignment start
            (0, 0i64, 1i64),  // invalid alignment start
            (0, 1i64, -1i64), // invalid alignment span
            (0, 1i64, 0i64),  // invalid alignment span
        ] {
            assert!(matches!(
                ReferenceSequenceContext::try_from(triplet),
                Err(e) if e.kind() == io::ErrorKind::InvalidData
            ));
        }

        Ok(())
    }
}
