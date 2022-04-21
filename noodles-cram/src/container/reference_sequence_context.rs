use std::cmp;

use noodles_core::Position;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceSequenceContext {
    Some(usize, Position, Position),
    None,
    Many,
}

impl ReferenceSequenceContext {
    pub fn update(
        &mut self,
        reference_sequence_id: Option<usize>,
        alignment_start: Option<Position>,
        alignment_end: Option<Position>,
    ) {
        *self = match (*self, reference_sequence_id, alignment_start, alignment_end) {
            (
                Self::Some(self_id, self_start, self_end),
                Some(record_id),
                Some(record_start),
                Some(record_end),
            ) => {
                if record_id == self_id {
                    let start = cmp::min(record_start, self_start);
                    let end = cmp::max(record_end, self_end);
                    Self::Some(record_id, start, end)
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

impl Default for ReferenceSequenceContext {
    fn default() -> Self {
        Self::None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_update() -> Result<(), noodles_core::position::TryFromIntError> {
        let mut context =
            ReferenceSequenceContext::Some(0, Position::try_from(8)?, Position::try_from(13)?);
        context.update(Some(0), Position::new(5), Position::new(21));
        assert_eq!(
            context,
            ReferenceSequenceContext::Some(0, Position::try_from(5)?, Position::try_from(21)?)
        );

        let mut context =
            ReferenceSequenceContext::Some(0, Position::try_from(8)?, Position::try_from(13)?);
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
}
