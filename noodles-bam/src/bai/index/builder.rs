use std::{cmp::Ordering, io, mem};

use noodles_csi::index::reference_sequence::bin::Chunk;
use noodles_sam::alignment::Record;

use super::{reference_sequence, Index};

/// A BAM index builder.
#[derive(Default)]
pub struct Builder {
    current_reference_sequence_id: usize,
    reference_sequence_builder: reference_sequence::Builder,
    reference_sequence_builders: Vec<reference_sequence::Builder>,
    unplaced_unmapped_record_count: u64,
}

impl Builder {
    /// Adds a record.
    ///
    /// The record must have an associated chunk denoting its start and end
    /// position in the file.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// use noodles_bgzf as bgzf;
    /// use noodles_csi::index::reference_sequence::bin::Chunk;
    /// use noodles_sam::alignment::Record;
    ///
    /// let mut builder = bai::Index::builder();
    ///
    /// let record = Record::default();
    /// let chunk = Chunk::new(
    ///     bgzf::VirtualPosition::from(233),
    ///     bgzf::VirtualPosition::from(377),
    /// );
    ///
    /// builder.add_record(&record, chunk);
    /// ```
    pub fn add_record(&mut self, record: &Record, chunk: Chunk) -> io::Result<()> {
        let (reference_sequence_id, start, end) = match (
            record.reference_sequence_id(),
            record.alignment_start(),
            record.alignment_end(),
        ) {
            (Some(reference_sequence_id), Some(start), Some(end)) => {
                (reference_sequence_id, start, end)
            }
            _ => {
                self.unplaced_unmapped_record_count += 1;
                return Ok(());
            }
        };

        match reference_sequence_id.cmp(&self.current_reference_sequence_id) {
            Ordering::Less => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                    "reference sequence ID ({}) appears after current reference sequence ID ({})",
                    reference_sequence_id, self.current_reference_sequence_id
                ),
                ))
            }
            Ordering::Equal => {}
            Ordering::Greater => self.add_reference_sequences_builders_until(reference_sequence_id),
        }

        self.reference_sequence_builder
            .add_record(start, end, record.flags(), chunk)
    }

    fn add_reference_sequences_builders_until(&mut self, reference_sequence_id: usize) {
        while self.current_reference_sequence_id < reference_sequence_id {
            let reference_sequence_builder = mem::take(&mut self.reference_sequence_builder);

            self.reference_sequence_builders
                .push(reference_sequence_builder);

            self.current_reference_sequence_id += 1;
        }
    }

    /// Builds a BAM index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// let index = bai::Index::builder().build(1);
    /// ```
    pub fn build(mut self, reference_sequence_count: usize) -> Index {
        if reference_sequence_count == 0 {
            return Index::new(Vec::new(), Some(self.unplaced_unmapped_record_count));
        }

        // SAFETY: `reference_sequence_count` is > 0.
        let last_reference_sequence_id = reference_sequence_count - 1;
        self.add_reference_sequences_builders_until(last_reference_sequence_id);

        self.reference_sequence_builders
            .push(self.reference_sequence_builder);

        let reference_sequences = self
            .reference_sequence_builders
            .into_iter()
            .map(|b| b.build())
            .collect();

        Index::new(
            reference_sequences,
            Some(self.unplaced_unmapped_record_count),
        )
    }
}

#[cfg(test)]
mod tests {
    use noodles_bgzf as bgzf;
    use noodles_core::Position;
    use noodles_csi::BinningIndex;
    use noodles_sam::record::Flags;

    use super::*;

    #[test]
    fn test_add_record_with_out_of_order_records() -> Result<(), Box<dyn std::error::Error>> {
        let mut builder = Builder::default();

        let record = Record::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(1)
            .set_alignment_start(Position::MIN)
            .set_cigar("4M".parse()?)
            .build();

        let chunk = Chunk::new(
            bgzf::VirtualPosition::from(55),
            bgzf::VirtualPosition::from(89),
        );

        builder.add_record(&record, chunk)?;

        let record = Record::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::MIN)
            .set_cigar("4M".parse()?)
            .build();

        let chunk = Chunk::new(
            bgzf::VirtualPosition::from(89),
            bgzf::VirtualPosition::from(144),
        );

        assert!(matches!(
            builder.add_record(&record, chunk),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        let mut builder = Builder::default();

        let record = Record::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(2)?)
            .set_cigar("4M".parse()?)
            .build();

        builder.add_record(
            &record,
            Chunk::new(
                bgzf::VirtualPosition::from(55),
                bgzf::VirtualPosition::from(89),
            ),
        )?;

        builder.add_record(
            &Record::default(),
            Chunk::new(
                bgzf::VirtualPosition::from(89),
                bgzf::VirtualPosition::from(144),
            ),
        )?;

        let index = builder.build(2);
        assert_eq!(index.reference_sequences().len(), 2);
        assert_eq!(index.unplaced_unmapped_record_count(), Some(1));

        Ok(())
    }
}
