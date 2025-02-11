use std::{io, mem, num};

use noodles_fasta as fasta;
use noodles_sam as sam;

use super::{slice, CompressionHeader, Container, Slice};
use crate::{io::writer::Options, Record};

const MAX_SLICE_COUNT: usize = 1;

#[derive(Debug)]
pub struct Builder {
    slice_builder: slice::Builder,
    slice_builders: Vec<slice::Builder>,
    record_counter: u64,
    base_count: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AddRecordError {
    InvalidRecordReadLength(num::TryFromIntError),
    ContainerFull(Record),
    SliceFull(Record),
}

impl Builder {
    pub fn new(record_counter: u64) -> Self {
        Self {
            slice_builder: Slice::builder(),
            slice_builders: Vec::new(),
            record_counter,
            base_count: 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.slice_builder.is_empty() && self.slice_builders.is_empty()
    }

    pub fn base_count(&self) -> u64 {
        self.base_count
    }

    #[allow(clippy::result_large_err)]
    pub fn add_record(&mut self, record: Record) -> Result<(), AddRecordError> {
        if self.slice_builders.len() >= MAX_SLICE_COUNT {
            return Err(AddRecordError::ContainerFull(record));
        }

        match self.slice_builder.add_record(record) {
            Ok(r) => {
                self.base_count += u64::try_from(r.read_length())
                    .map_err(AddRecordError::InvalidRecordReadLength)?;
                Ok(())
            }
            Err(e) => match e {
                slice::builder::AddRecordError::SliceFull(r) => {
                    let slice_builder = mem::take(&mut self.slice_builder);
                    self.slice_builders.push(slice_builder);
                    Err(AddRecordError::SliceFull(r))
                }
            },
        }
    }

    pub fn build(
        mut self,
        options: &Options,
        reference_sequence_repository: &fasta::Repository,
        header: &sam::Header,
    ) -> io::Result<Container> {
        if !self.slice_builder.is_empty() {
            self.slice_builders.push(self.slice_builder);
        }

        let mut options = options.clone();

        if self
            .slice_builders
            .iter()
            .any(|b| b.reference_sequence_context().is_many())
        {
            options.encode_alignment_start_positions_as_deltas = false;
        }

        let compression_header = build_compression_header(&options, &self.slice_builders);

        let record_counter = self.record_counter;
        let slices = self
            .slice_builders
            .into_iter()
            .map(|builder| {
                builder.build(
                    &options.block_content_encoder_map,
                    reference_sequence_repository,
                    header,
                    &compression_header,
                    record_counter,
                )
            })
            .collect::<Result<_, _>>()?;

        Ok(Container {
            compression_header,
            slices,
        })
    }
}

fn build_compression_header(
    options: &Options,
    slice_builders: &[slice::Builder],
) -> CompressionHeader {
    let mut compression_header_builder = CompressionHeader::builder();
    compression_header_builder.apply_options(options);

    for slice_builder in slice_builders {
        for record in slice_builder.records() {
            compression_header_builder.update(record);
        }
    }

    compression_header_builder.build()
}
