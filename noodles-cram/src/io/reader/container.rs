mod block;
pub mod compression_header;
pub mod header;
pub mod slice;

use std::{
    io::{self, Read},
    iter,
};

use bytes::BytesMut;

pub(crate) use self::block::Block;
use self::header::read_header;
pub use self::{block::read_block, compression_header::get_compression_header, slice::read_slice};
use crate::container::{CompressionHeader, Header, Slice};

/// A CRAM container.
#[derive(Default)]
pub struct Container {
    pub(crate) header: Header,
    pub(crate) src: BytesMut,
}

impl Container {
    /// Returns the header.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Returns the compression header.
    pub fn compression_header(&mut self) -> io::Result<CompressionHeader> {
        let end = self
            .header
            .landmarks
            .first()
            .copied()
            .unwrap_or(self.src.len());

        let src = self.src.split_to(end);
        let mut buf = src.freeze();

        get_compression_header(&mut buf)
    }

    /// Returns the iterator over slices.
    pub fn slices(&mut self) -> impl Iterator<Item = io::Result<Slice>> + '_ {
        let landmarks = self.header.landmarks.clone();
        let mut i = 0;

        iter::from_fn(move || {
            if i < landmarks.len() - 1 {
                let end = landmarks[i + 1];
                i += 1;
                let src = self.src.split_to(end);
                let mut buf = src.freeze();
                Some(read_slice(&mut buf))
            } else if i < landmarks.len() {
                let src = self.src.split();
                let mut buf = src.freeze();
                Some(read_slice(&mut buf))
            } else {
                None
            }
        })
    }
}

pub fn read_container<R>(reader: &mut R, container: &mut Container) -> io::Result<usize>
where
    R: Read,
{
    match read_header(reader, &mut container.header)? {
        0 => Ok(0),
        len => {
            container.src.resize(len, 0);
            reader.read_exact(&mut container.src)?;
            Ok(len)
        }
    }
}
