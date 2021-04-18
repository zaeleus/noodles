//! A pileup engine over sets of reads.
// TODOs and want list
// - [] add the query base to the `RecordAlignment`
// - [] add more constructors for Pileup
// - [] Either make sure a Reader allows for filtering, or allow for filtering callback filter reads before they enter the pileup
// - [] Add a nice debug repr of a PileupColumn
// - [] Add a method to FixMates as in perbase / make sure enough info is present to do this
// - [] Add tests
// - [] Add examples to docs
//
// Other information
// htslib pileup: https://docs.rs/rust-htslib/0.36.0/rust_htslib/bam/pileup/struct.Pileup.html
// htslib alignment: https://docs.rs/rust-htslib/0.36.0/rust_htslib/bam/pileup/struct.Alignment.html

use std::{
    collections::VecDeque,
    convert::TryFrom,
    io::{self, ErrorKind, Read, Seek},
    iter::Peekable,
};

use super::{reader::Query, record::cigar::Op, Record};
use noodles_sam::record::cigar::op::Kind;

/// Wrapper of a CIGAR Op to track where in the Op we are
#[derive(Debug)]
struct OpWrapper {
    /// The Cigar Op.
    op: Op,
    /// The remaining number of the given operation. This will be decremented as operations are consumed.
    len: u32,
}
/// [`AugmentedRecord`] wraps a record and pre-computes several fields based on the underlying record.
/// It also holds onto the state of where in the read we are.
#[derive(Debug)]
struct AugmentedRecord {
    /// The wrapped underlying [`Record`].
    record: Record,
    /// The reads start position on the reference.
    start: i32,
    /// The reads end position on the reference.
    ref_end: i32,
    /// A queue of wrapped Cigar Ops. This will be consumed as the read is iterated over.
    cigar_ops: VecDeque<OpWrapper>,
}

impl AugmentedRecord {
    /// Get the current cigar op.
    fn consume_op(&mut self) -> Op {
        let mut front_op = &mut self.cigar_ops[0];
        front_op.len -= 1;
        let op = front_op.op;
        // If we've exhausted this op, pop from queue
        if front_op.len == 0 {
            let _ = self.cigar_ops.pop_front();
        }
        op
    }

    /// Skip any clipping or padding at the start of the read.
    fn skip_clipping(&mut self) {
        let mut num_pop = 0;
        let mut iter = self.cigar_ops.iter().map(|c_op| c_op.op.kind());
        while let Some(Kind::SoftClip | Kind::HardClip | Kind::Pad) = iter.next() {
            num_pop += 1;
        }
        for _ in 0..num_pop {
            let _ = self.cigar_ops.pop_front();
        }
    }

    /// Create a [`RecordAlignment`] based on the current position.
    fn as_record_alignment(&mut self) -> RecordAlignment {
        RecordAlignment {
            op: self.consume_op(),
        }
    }
}

impl<'a> TryFrom<Record> for AugmentedRecord {
    type Error = io::Error;

    /// Create an [`AugmentedRecord`] from a [`Record`].
    fn try_from(record: Record) -> Result<Self, Self::Error> {
        let cigar = record.cigar();
        let mut cigar_ops = VecDeque::new();
        for op in cigar.ops() {
            let op = op?;
            cigar_ops.push_back(OpWrapper { op, len: op.len() })
        }
        // NB: unwrap is safe since only mapped reads are present here
        let start = record.position().unwrap().into();
        let ref_end = start + cigar.reference_len()? as i32;
        Ok(AugmentedRecord {
            record,
            start,
            ref_end,
            cigar_ops,
        })
    }
}

/// The information about a position in a specific Record
#[derive(Debug)]
pub struct RecordAlignment {
    /// The cigar [`Op`] for this read at this position
    pub op: Op,
}

/// The pileup over a specific position.
#[derive(Debug)]
pub struct PileupColumn {
    /// 1-based position in the reference
    pub pos: i32,
    /// The total depth at this position.
    pub depth: u32,
    /// Number of insertions that start to the right of this position. (does not count toward depth).
    pub ins: u32,
    /// Number of deletions at this position (counts towards depth).
    pub del: u32,
    /// Number of refskips at this position. Not counted toward depth.
    /// Note: this differs from how htslib pileup computes depths
    pub ref_skip: u32,
    /// The alignments at this position.
    pub alignments: Vec<RecordAlignment>,
}

impl PileupColumn {
    /// Build a Pileup over `pos` from an iterator over [`RecordAlignment`]s.
    fn from_alignments(alignments: impl Iterator<Item = RecordAlignment>, pos: i32) -> Self {
        let mut depth = 0;
        let mut ins = 0;
        let mut del = 0;
        let mut ref_skip = 0;
        let mut alignment_records = Vec::with_capacity(alignments.size_hint().0);
        for alignment in alignments {
            match alignment.op.kind() {
                Kind::Match | Kind::SeqMatch | Kind::SeqMismatch => {
                    depth += 1;
                }
                Kind::Insertion => {
                    ins += 1;
                }
                Kind::Deletion => {
                    del += 1;
                    depth += 1;
                }
                Kind::Skip => {
                    ref_skip += 1;
                }
                // Should never see?
                Kind::SoftClip | Kind::HardClip | Kind::Pad => unreachable!(),
            }
            alignment_records.push(alignment)
        }
        PileupColumn {
            pos,
            depth,
            ins,
            del,
            ref_skip,
            alignments: alignment_records,
        }
    }
}

/// Hold state information to generate column pileups over a query regions
pub struct Pileup<'a> {
    /// A [`Peekable`] iterator of coordinate sorted [`AugmentedRecord`]s.
    iterator: Peekable<Box<dyn Iterator<Item = Result<AugmentedRecord, io::Error>> + 'a>>,
    /// The queue of reads that overlap the current position.
    queue: VecDeque<AugmentedRecord>,
    /// The current reference position.
    curr_position: i32,
}

impl<'a> Pileup<'a> {
    /// Create a new [`Pileup`] from a [`Query`]
    pub fn new<R: Read + Seek>(query: Query<'a, R>) -> Result<Self, io::Error> {
        // Wrap the query iterator in a transformer to create AugmentedRecords
        let mut reads = (Box::new(query.map(|read| {
            let read = read.unwrap();
            AugmentedRecord::try_from(read).map(|mut rec| {
                rec.skip_clipping();
                rec
            })
        }))
            as Box<dyn Iterator<Item = Result<AugmentedRecord, io::Error>>>)
            .peekable();

        // Try to get the first position in read set
        let curr_pos = if let Some(read) = reads.peek() {
            let read = read
                .as_ref()
                .map_err(|e| io::Error::new(ErrorKind::Other, e.to_string()))?;
            read.start
        } else {
            1
        };

        Ok(Self {
            iterator: reads,
            curr_position: curr_pos,
            queue: VecDeque::new(),
        })
    }

    /// Fill the queue with reads until the reads no longer overlap current position
    fn fill_queue(&mut self) -> Result<(), io::Error> {
        while let Some(read) = self.iterator.peek() {
            let read = read
                .as_ref()
                .map_err(|e| io::Error::new(ErrorKind::Other, e.to_string()))?;

            // NB: below this point we know we can unwrap since we successfully peeked ahead
            // NB: Skip unmapped reads
            if read.record.flags().is_unmapped() {
                let _ = self.iterator.next().unwrap();
                continue;
            }
            if read.start <= self.curr_position && read.ref_end > self.curr_position {
                self.queue.push_back(self.iterator.next().unwrap()?);
            } else if self.queue.is_empty() {
                // if the queue is empty and the next read exists but doesn't overlap curr_position, jump ahead
                self.curr_position = read.start;
                self.queue.push_back(self.iterator.next().unwrap()?);
            } else {
                break;
            }
        }
        Ok(())
    }

    /// Create a [`PileupColumn`] based on the current queue and current position.
    fn build_column(&mut self) -> PileupColumn {
        let alignments = self.queue.iter_mut().map(|rec| rec.as_record_alignment());
        PileupColumn::from_alignments(alignments, self.curr_position)
    }

    /// Remove reads that no longer overlap the current position.
    fn evict_queue(&mut self) {
        let curr_pos = self.curr_position;

        // TODO: make queue a vec? / do what sambamba does and overwrite values from start?
        for i in (0..self.queue.len()).rev() {
            let read = &self.queue[i];
            if !(read.start <= curr_pos && read.ref_end > curr_pos) {
                self.queue.swap_remove_back(i);
            }
        }
    }
}

impl<'a> Iterator for Pileup<'a> {
    type Item = Result<PileupColumn, io::Error>;
    /// Get the next positions [`PileupColumn`].
    fn next(&mut self) -> Option<Self::Item> {
        if let Err(err) = self.fill_queue() {
            return Some(Err(err));
        }
        if !self.queue.is_empty() {
            let pileup_column = self.build_column();
            self.curr_position += 1;
            self.evict_queue();
            Some(Ok(pileup_column))
        } else {
            None
        }
    }
}
