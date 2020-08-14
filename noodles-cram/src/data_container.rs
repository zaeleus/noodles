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
    pub fn builder() -> Builder {
        Builder::default()
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

        let compression_header =
            CompressionHeader::try_from(&blocks[0]).expect("missing compression header");

        let mut start = 1;

        let slices_len = container.header().landmarks().len();
        let mut slices = Vec::with_capacity(slices_len);

        for _ in 0..slices_len {
            let slice = Slice::try_from(&blocks[start..]).expect("missing slice");
            let block_count = slice.header().block_count();

            slices.push(slice);

            start += block_count as usize;
        }

        Ok(DataContainer {
            compression_header,
            slices,
        })
    }
}
