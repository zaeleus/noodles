pub mod builder;

pub use self::builder::Builder;

use std::{convert::TryFrom, io};

use super::{
    container::{CompressionHeader, Slice},
    Container,
};

pub struct DataContainer {
    compression_header: CompressionHeader,
    slices: Vec<Slice>,
}

impl DataContainer {
    pub fn builder(record_counter: i64) -> Builder {
        Builder::new(record_counter)
    }

    pub fn compression_header(&self) -> &CompressionHeader {
        &self.compression_header
    }

    pub fn slices(&self) -> &[Slice] {
        &self.slices
    }
}

impl TryFrom<Container> for DataContainer {
    type Error = io::Error;

    fn try_from(container: Container) -> Result<Self, Self::Error> {
        let blocks = container.blocks();

        let compression_header = blocks
            .first()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "missing compression header block",
                )
            })
            .and_then(CompressionHeader::try_from)?;

        let mut start = 1;

        let slices_len = container.header().landmarks().len();
        let mut slices = Vec::with_capacity(slices_len);

        for _ in 0..slices_len {
            let slice = Slice::try_from(&blocks[start..])?;

            // (core data block + external blocks) + header block
            let block_count = slice.header().block_count() + 1;

            slices.push(slice);

            start += block_count as usize;
        }

        Ok(DataContainer {
            compression_header,
            slices,
        })
    }
}
