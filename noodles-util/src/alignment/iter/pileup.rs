use std::{collections::VecDeque, io};

use noodles_core::Position;
use noodles_sam::{
    self as sam, Header,
    alignment::{Record, record::Flags},
};

type ActiveWindowRange = (Position, Position);

#[derive(Debug)]
enum State {
    Empty,
    Pile(ActiveWindowRange),
    Pop(ActiveWindowRange),
    Drain,
    Done,
}

/// A pileup iterator.
///
/// This takes an iterator of coordinate-sorted records and emits reference sequence column
/// statistics.
pub struct Pileup<'h, I> {
    header: &'h Header,
    records: I,
    state: State,
    position: Position,
    window: VecDeque<u64>,
    next_record: Option<Box<dyn Record>>,
}

impl<'h, I> Pileup<'h, I>
where
    I: Iterator<Item = io::Result<Box<dyn Record>>>,
{
    /// Creates a pileup iterator.
    ///
    /// The given iterator must be coordinate-sorted on a single reference sequence.
    pub fn new(header: &'h Header, records: I) -> Self {
        Self {
            header,
            records,
            state: State::Empty,
            position: Position::MIN,
            window: VecDeque::new(),
            next_record: None,
        }
    }

    fn initialize(&mut self) -> io::Result<Option<ActiveWindowRange>> {
        if self.next_record.is_none() {
            for result in &mut self.records {
                let record = result?;
                let flags = record.flags()?;

                if filter(flags) {
                    continue;
                }

                self.next_record = Some(record);

                break;
            }
        }

        if let Some(record) = self.next_record.take() {
            let (_, start, end) = alignment_context(self.header, &record)?;
            self.position = start;
            pile_record(&mut self.window, start, end, &record)?;
            Ok(Some((start, end)))
        } else {
            Ok(None)
        }
    }

    fn pile_records(
        &mut self,
        active_window_range: ActiveWindowRange,
    ) -> io::Result<Option<ActiveWindowRange>> {
        let (mut active_window_start, mut active_window_end) = active_window_range;

        if let Some(record) = self.next_record.take() {
            let (_, start, end) = alignment_context(self.header, &record)?;
            pile_record(&mut self.window, start, end, &record)?;
            active_window_end = end.max(active_window_end);
        }

        while let Some(record) = self.records.next().transpose()? {
            let flags = record.flags()?;

            if filter(flags) {
                continue;
            }

            let (_, start, end) = alignment_context(self.header, &record)?;

            if start > active_window_end {
                self.next_record = Some(record);
                return Ok(None);
            } else if start > active_window_start {
                self.next_record = Some(record);
                active_window_start = start;
                return Ok(Some((active_window_start, active_window_end)));
            }

            pile_record(&mut self.window, start, end, &record)?;
            active_window_end = end.max(active_window_end);
        }

        Ok(None)
    }

    fn pop_front_full(&mut self) -> Option<(Position, u64)> {
        let position = self.position;
        let record = self.window.pop_front()?;

        self.position = self
            .position
            .checked_add(1)
            .expect("attempt to add with overflow");

        Some((position, record))
    }
}

impl<I> Iterator for Pileup<'_, I>
where
    I: Iterator<Item = io::Result<Box<dyn Record>>>,
{
    type Item = io::Result<(Position, u64)>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.state = match self.state {
                State::Empty => match self.initialize() {
                    Ok(None) => State::Done,
                    Ok(Some(active_window_range)) => State::Pile(active_window_range),
                    Err(e) => return Some(Err(e)),
                },
                State::Pile(active_window_range) => match self.pile_records(active_window_range) {
                    Ok(None) => State::Drain,
                    Ok(Some(next_active_window_range)) => State::Pop(next_active_window_range),
                    Err(e) => return Some(Err(e)),
                },
                State::Pop((active_window_start, active_window_end)) => {
                    if self.position < active_window_start {
                        // SAFETY: active_window_start - self.position < self.window.len()
                        let value = self.pop_front_full().unwrap();
                        return Some(Ok(value));
                    } else {
                        State::Pile((active_window_start, active_window_end))
                    }
                }
                State::Drain => match self.pop_front_full() {
                    Some(value) => return Some(Ok(value)),
                    None => State::Empty,
                },
                State::Done => return None,
            }
        }
    }
}

fn alignment_context<R>(header: &Header, record: &R) -> io::Result<(usize, Position, Position)>
where
    R: Record,
{
    match (
        record.reference_sequence_id(header).transpose()?,
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ) {
        (Some(id), Some(start), Some(end)) => Ok((id, start, end)),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing reference sequence ID or alignment start",
        )),
    }
}

fn filter(flags: Flags) -> bool {
    flags.is_unmapped() || flags.is_secondary() || flags.is_qc_fail() || flags.is_duplicate()
}

fn pile_record<R>(
    window: &mut VecDeque<u64>,
    start: Position,
    end: Position,
    record: &R,
) -> io::Result<()>
where
    R: Record,
{
    let span = usize::from(end) - usize::from(start) + 1;

    if span > window.len() {
        window.resize(span, 0);
    }

    let cigar = record.cigar();
    pile(window, start, start, &cigar)
}

fn pile<C>(
    window: &mut VecDeque<u64>,
    offset: Position,
    start: Position,
    cigar: &C,
) -> io::Result<()>
where
    C: sam::alignment::record::Cigar,
{
    use sam::alignment::record::cigar::op::Kind;

    let offset = usize::from(offset) - 1;
    let start = usize::from(start) - 1;
    let mut i = start - offset;

    for result in cigar.iter() {
        let op = result?;

        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let end = i + op.len();

                for depth in window.range_mut(i..end) {
                    *depth += 1;
                }

                i = end;
            }
            Kind::Deletion | Kind::Skip => i += op.len(),
            _ => {}
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use super::*;
    use sam::alignment::RecordBuf;

    #[test]
    fn test_next() -> Result<(), Box<dyn std::error::Error>> {
        use sam::{
            alignment::record::cigar::{Op, op::Kind},
            header::record::value::{Map, map::ReferenceSequence},
        };

        // 1 2 3 4 5 6 7 8 9
        //   [   ]
        //   [     ]
        //     [ ]
        //       [ ]
        //             [ ]
        //             [   ]
        let records: Vec<_> = [
            (
                0,
                Position::try_from(2)?,
                [Op::new(Kind::Match, 3)].into_iter().collect(),
            ),
            (
                0,
                Position::try_from(2)?,
                [Op::new(Kind::Match, 4)].into_iter().collect(),
            ),
            (
                0,
                Position::try_from(3)?,
                [Op::new(Kind::Match, 2)].into_iter().collect(),
            ),
            (
                0,
                Position::try_from(4)?,
                [Op::new(Kind::Match, 2)].into_iter().collect(),
            ),
            (
                0,
                Position::try_from(7)?,
                [Op::new(Kind::Match, 2)].into_iter().collect(),
            ),
            (
                0,
                Position::try_from(7)?,
                [Op::new(Kind::Match, 3)].into_iter().collect(),
            ),
        ]
        .into_iter()
        .map(|(reference_sequence_id, position, cigar)| {
            RecordBuf::builder()
                .set_flags(Flags::empty())
                .set_reference_sequence_id(reference_sequence_id)
                .set_alignment_start(position)
                .set_cigar(cigar)
                .build()
        })
        .map(|record| Ok(Box::new(record) as Box<dyn Record>))
        .collect();

        let header = Header::builder()
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(NonZero::<usize>::MAX))
            .build();

        let pileup = Pileup::new(&header, records.into_iter());
        let actual: Vec<_> = pileup.collect::<Result<_, _>>()?;

        let expected = [
            (Position::try_from(2)?, 2),
            (Position::try_from(3)?, 3),
            (Position::try_from(4)?, 4),
            (Position::try_from(5)?, 2),
            (Position::try_from(7)?, 2),
            (Position::try_from(8)?, 2),
            (Position::try_from(9)?, 1),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
