use std::{collections::VecDeque, io};

use noodles_core::Position;

use crate::{alignment::Record, record::Flags, Header};

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
                let flags = try_to_flags(record.flags().as_ref())?;

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
            pile_record(&mut self.window, start, end, self.header, &record)?;
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
            pile_record(&mut self.window, start, end, self.header, &record)?;
            active_window_end = end.max(active_window_end);
        }

        while let Some(record) = self.records.next().transpose()? {
            let flags = try_to_flags(record.flags().as_ref())?;

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

            pile_record(&mut self.window, start, end, self.header, &record)?;
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

impl<'a, I> Iterator for Pileup<'a, I>
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

fn try_to_flags(flags: &dyn crate::alignment::record::Flags) -> io::Result<Flags> {
    flags.try_to_u16().map(Flags::from)
}

fn try_to_reference_sequence_id<I>(reference_sequence_id: I) -> io::Result<usize>
where
    I: crate::alignment::record::ReferenceSequenceId,
{
    reference_sequence_id.try_to_usize()
}

fn try_to_position<P>(position: P) -> io::Result<Position>
where
    P: crate::alignment::record::Position,
{
    Position::try_from(&position as &dyn crate::alignment::record::Position)
}

fn alignment_context<R>(header: &Header, record: &R) -> io::Result<(usize, Position, Position)>
where
    R: Record,
{
    match (
        record.reference_sequence_id(header),
        record.alignment_start(),
        record.alignment_end(header),
    ) {
        (Some(id), Some(start), Some(end)) => {
            let id = try_to_reference_sequence_id(id)?;
            let start = try_to_position(start)?;
            let end = end?;
            Ok((id, start, end))
        }
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
    header: &Header,
    record: &R,
) -> io::Result<()>
where
    R: Record,
{
    let span = usize::from(end) - usize::from(start) + 1;

    if span > window.len() {
        window.resize(span, 0);
    }

    let cigar = record.cigar(header);
    pile(window, start, start, &cigar)
}

fn pile<C>(
    window: &mut VecDeque<u64>,
    offset: Position,
    start: Position,
    cigar: &C,
) -> io::Result<()>
where
    C: crate::alignment::record::Cigar,
{
    use crate::record::cigar::op::Kind;

    let offset = usize::from(offset) - 1;
    let start = usize::from(start) - 1;
    let mut i = start - offset;

    for result in cigar.iter() {
        let (kind, len) = result?;

        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let end = i + len;

                for depth in window.range_mut(i..end) {
                    *depth += 1;
                }

                i = end;
            }
            Kind::Deletion | Kind::Skip => i += len,
            _ => {}
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use super::*;
    use crate::alignment::RecordBuf;

    #[test]
    fn test_next() -> Result<(), Box<dyn std::error::Error>> {
        use crate::header::record::value::{map::ReferenceSequence, Map};

        // 1 2 3 4 5 6 7 8 9
        //   [   ]
        //   [     ]
        //     [ ]
        //       [ ]
        //             [ ]
        //             [   ]
        let records: Vec<_> = [
            (0, Position::try_from(2)?, "3M".parse()?),
            (0, Position::try_from(2)?, "4M".parse()?),
            (0, Position::try_from(3)?, "2M".parse()?),
            (0, Position::try_from(4)?, "2M".parse()?),
            (0, Position::try_from(7)?, "2M".parse()?),
            (0, Position::try_from(7)?, "3M".parse()?),
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
            .add_reference_sequence(
                "sq0".parse()?,
                Map::<ReferenceSequence>::new(NonZeroUsize::MAX),
            )
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
