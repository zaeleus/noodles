use std::slice;

use noodles_core::Position;
use noodles_sam::alignment::record::cigar::op::Kind;

use crate::record::Feature;

/// An iterator over features as CIGAR operations.
pub struct Cigar<'a> {
    features: slice::Iter<'a, Feature>,
    read_length: usize,
    read_position: Position,
    next_op: Option<(Kind, usize)>,
}

impl<'a> Cigar<'a> {
    pub(super) fn new(features: &'a [Feature], read_length: usize) -> Self {
        Self {
            features: features.iter(),
            read_length,
            read_position: Position::MIN,
            next_op: None,
        }
    }

    fn consume_read(&mut self, len: usize) {
        self.read_position = self
            .read_position
            .checked_add(len)
            .expect("attempt to add with overflow");
    }
}

impl<'a> Iterator for Cigar<'a> {
    type Item = (Kind, usize);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(op) = self.next_op.take() {
            return Some(op);
        }

        let Some(feature) = self.features.next() else {
            if usize::from(self.read_position) <= self.read_length {
                let len = self.read_length - usize::from(self.read_position) + 1;
                self.consume_read(len);
                return Some((Kind::Match, len));
            } else {
                return None;
            }
        };

        if feature.position() > self.read_position {
            let len = usize::from(feature.position()) - usize::from(self.read_position);
            self.read_position = feature.position();
            self.next_op = Some((Kind::Match, len));
        }

        let (kind, len) = match feature {
            Feature::Substitution(..) => (Kind::Match, 1),
            Feature::Insertion(_, bases) => (Kind::Insertion, bases.len()),
            Feature::Deletion(_, len) => (Kind::Deletion, *len),
            Feature::InsertBase(..) => (Kind::Insertion, 1),
            Feature::ReferenceSkip(_, len) => (Kind::Skip, *len),
            Feature::SoftClip(_, bases) => (Kind::SoftClip, bases.len()),
            Feature::Padding(_, len) => (Kind::Pad, *len),
            Feature::HardClip(_, len) => (Kind::HardClip, *len),
            _ => todo!(),
        };

        if kind.consumes_read() {
            self.consume_read(len);
        }

        match self.next_op.replace((kind, len)) {
            Some(op) => Some(op),
            None => self.next_op.take(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::{feature::substitution, Features};

    #[test]
    fn test_next() -> Result<(), noodles_core::position::TryFromIntError> {
        fn t(features: &Features, read_length: usize, expected: &[(Kind, usize)]) {
            let cigar = Cigar::new(features, read_length);
            let actual: Vec<_> = cigar.collect();
            assert_eq!(actual, expected);
        }

        let features = Features::default();
        t(&features, 4, &[(Kind::Match, 4)]);

        let features = Features::from(vec![Feature::SoftClip(
            Position::try_from(1)?,
            vec![b'A', b'T'],
        )]);
        t(&features, 4, &[(Kind::SoftClip, 2), (Kind::Match, 2)]);

        let features = Features::from(vec![Feature::SoftClip(Position::try_from(4)?, vec![b'G'])]);
        t(&features, 4, &[(Kind::Match, 3), (Kind::SoftClip, 1)]);

        let features = Features::from(vec![Feature::HardClip(Position::try_from(1)?, 2)]);
        t(&features, 4, &[(Kind::HardClip, 2), (Kind::Match, 4)]);

        // FIXME
        let features = Features::from(vec![
            Feature::SoftClip(Position::try_from(1)?, vec![b'A']),
            Feature::Substitution(Position::try_from(3)?, substitution::Value::Code(0)),
        ]);
        t(
            &features,
            4,
            &[
                (Kind::SoftClip, 1),
                (Kind::Match, 1),
                (Kind::Match, 1),
                (Kind::Match, 1),
            ],
        );

        // FIXME
        let features = Features::from(vec![Feature::Substitution(
            Position::try_from(2)?,
            substitution::Value::Code(0),
        )]);
        t(
            &features,
            4,
            &[(Kind::Match, 1), (Kind::Match, 1), (Kind::Match, 2)],
        );

        Ok(())
    }
}
