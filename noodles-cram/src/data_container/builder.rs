use std::{io, mem};

use noodles_fasta as fasta;

use crate::{
    container::{
        compression_header,
        slice::{self, Slice},
        CompressionHeader,
    },
    Record,
};

use super::DataContainer;

#[derive(Debug)]
pub struct Builder {
    compression_header_builder: compression_header::Builder,
    slice_builder: slice::Builder,
    slice_builders: Vec<slice::Builder>,
    record_counter: i64,
    base_count: i64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AddRecordError {
    ContainerFull(Record),
}

impl Builder {
    pub fn new(record_counter: i64) -> Self {
        Self {
            compression_header_builder: CompressionHeader::builder(),
            slice_builder: Slice::builder(),
            slice_builders: Vec::new(),
            record_counter,
            base_count: 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.slice_builder.is_empty() && self.slice_builders.is_empty()
    }

    pub fn base_count(&self) -> i64 {
        self.base_count
    }

    pub fn add_record(
        &mut self,
        reference_sequence: &[u8],
        record: Record,
    ) -> Result<(), AddRecordError> {
        match self.slice_builder.add_record(record) {
            Ok(r) => {
                self.compression_header_builder
                    .update(reference_sequence, r);

                self.base_count += i64::from(r.read_length());

                Ok(())
            }
            Err(e) => match e {
                slice::builder::AddRecordError::ReferenceSequenceIdMismatch(r) => {
                    let slice_builder = mem::replace(&mut self.slice_builder, Slice::builder());
                    self.slice_builders.push(slice_builder);
                    Err(AddRecordError::ContainerFull(r))
                }
            },
        }
    }

    pub fn build(mut self, reference_sequences: &[fasta::Record]) -> io::Result<DataContainer> {
        if !self.slice_builder.is_empty() {
            self.slice_builders.push(self.slice_builder);
        }

        let compression_header = self.compression_header_builder.build();

        let record_counter = self.record_counter;
        let slices = self
            .slice_builders
            .into_iter()
            .map(|builder| builder.build(reference_sequences, &compression_header, record_counter))
            .collect::<Result<_, _>>()?;

        Ok(DataContainer {
            compression_header,
            slices,
        })
    }
}
