use std::{io, mem};

use noodles_fasta as fasta;

use crate::{
    container::{
        compression_header,
        slice::{self, Slice},
    },
    Record,
};

use super::DataContainer;

#[derive(Debug, Default)]
pub struct Builder {
    compression_header_builder: compression_header::Builder,
    slice_builder: slice::Builder,
    slice_builders: Vec<slice::Builder>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AddRecordError {
    ContainerFull(Record),
}

impl Builder {
    pub fn add_record(
        &mut self,
        reference_sequence: &[u8],
        record: Record,
    ) -> Result<(), AddRecordError> {
        match self.slice_builder.add_record(record) {
            Ok(r) => {
                self.compression_header_builder
                    .update(reference_sequence, r);
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

        let slices = self
            .slice_builders
            .into_iter()
            .map(|builder| builder.build(reference_sequences, &compression_header))
            .collect::<Result<_, _>>()?;

        Ok(DataContainer {
            compression_header,
            slices,
        })
    }
}
