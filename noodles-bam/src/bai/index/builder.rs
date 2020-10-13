use std::io;

use crate::{record::ReferenceSequenceId, Record};

use super::{
    reference_sequence::{self, bin::Chunk},
    Index, ReferenceSequence,
};

/// A BAM index builder.
#[derive(Default)]
pub struct Builder {
    current_reference_sequence_id: ReferenceSequenceId,
    reference_sequences_builders: Vec<reference_sequence::Builder>,
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
    /// use noodles_bam::{self as bam, bai::{self, index::reference_sequence::bin::Chunk}};
    /// use noodles_bgzf as bgzf;
    ///
    /// let mut builder = bai::Index::builder();
    ///
    /// let record = bam::Record::default();
    /// let chunk = Chunk::new(
    ///     bgzf::VirtualPosition::from(233),
    ///     bgzf::VirtualPosition::from(377),
    /// );
    ///
    /// builder.add_record(&record, chunk);
    /// ```
    pub fn add_record(&mut self, record: &Record, chunk: Chunk) -> io::Result<()> {
        if record.position().is_none() {
            self.unplaced_unmapped_record_count += 1;
            return Ok(());
        }

        if record.reference_sequence_id() != self.current_reference_sequence_id {
            self.add_reference_sequences_builders_until(record.reference_sequence_id());
        }

        let reference_sequence_builder = self.reference_sequences_builders.last_mut().unwrap();

        reference_sequence_builder.add_record(record, chunk)
    }

    fn add_reference_sequences_builders_until(
        &mut self,
        reference_sequence_id: ReferenceSequenceId,
    ) {
        let id = i32::from(reference_sequence_id);
        let mut current_id = i32::from(self.current_reference_sequence_id);

        while current_id < id {
            self.reference_sequences_builders
                .push(ReferenceSequence::builder());
            current_id += 1;
        }

        self.current_reference_sequence_id = reference_sequence_id;
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
        let last_reference_sequence_id =
            ReferenceSequenceId::from((reference_sequence_count - 1) as i32);
        self.add_reference_sequences_builders_until(last_reference_sequence_id);

        let reference_sequences = self
            .reference_sequences_builders
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
    use std::convert::TryFrom;

    use noodles_bgzf as bgzf;
    use noodles_sam::{
        self as sam,
        header::ReferenceSequences,
        record::{Flags, Position},
    };

    use super::*;

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        let mut reference_sequences = ReferenceSequences::default();
        reference_sequences.insert(
            String::from("sq0"),
            sam::header::ReferenceSequence::new(String::from("sq0"), 8),
        );
        reference_sequences.insert(
            String::from("sq1"),
            sam::header::ReferenceSequence::new(String::from("sq1"), 13),
        );

        let mut builder = Builder::default();

        let record = Record::try_from_sam_record(
            &reference_sequences,
            &sam::Record::builder()
                .set_flags(Flags::empty())
                .set_reference_sequence_name("sq0".parse()?)
                .set_position(Position::try_from(2)?)
                .set_cigar("4M".parse()?)
                .build(),
        )?;

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

        let index = builder.build(reference_sequences.len());
        assert_eq!(index.reference_sequences().len(), 2);
        assert_eq!(index.unplaced_unmapped_read_count(), Some(1));

        Ok(())
    }
}
