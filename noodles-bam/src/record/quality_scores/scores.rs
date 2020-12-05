use std::convert::TryFrom;

use noodles_sam::record::quality_scores::{score, Score};

/// An iterator over quality scores.
///
/// This is created by calling [`super::QualityScores::chars`].
pub struct Scores<I> {
    iter: I,
}

impl<I> Scores<I> {
    pub(crate) fn new(iter: I) -> Self {
        Self { iter }
    }
}

impl<'a, I> Iterator for Scores<I>
where
    I: Iterator<Item = &'a u8>,
{
    type Item = Result<Score, score::TryFromUByteError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().cloned().map(Score::try_from)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}
