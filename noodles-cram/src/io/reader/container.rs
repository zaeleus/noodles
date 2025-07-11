pub(crate) mod block;
pub mod compression_header;
pub mod header;
pub mod slice;

use std::{
    io::{self, Read},
    iter,
};

use self::{block::read_block_as, header::read_header};
pub use self::{compression_header::read_compression_header, slice::Slice, slice::read_slice};
use crate::container::{CompressionHeader, Header};

/// A CRAM container.
#[derive(Default)]
pub struct Container {
    pub(crate) header: Header,
    pub(crate) src: Vec<u8>,
}

impl Container {
    /// Returns the header.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Returns the compression header.
    pub fn compression_header(&self) -> io::Result<CompressionHeader> {
        let end = self
            .header
            .landmarks
            .first()
            .copied()
            .unwrap_or(self.src.len());

        let mut src = &self.src[..end];

        read_compression_header(&mut src)
    }

    /// Returns the iterator over slices.
    pub fn slices(&self) -> impl Iterator<Item = io::Result<Slice<'_>>> + '_ {
        let landmarks = &self.header.landmarks;
        let mut i = 0;

        iter::from_fn(move || {
            if i < landmarks.len() - 1 {
                let (start, end) = (landmarks[i], landmarks[i + 1]);
                i += 1;
                let mut src = &self.src[start..end];
                Some(read_slice(&mut src))
            } else if i < landmarks.len() {
                let start = landmarks[i];
                i += 1;
                let mut src = &self.src[start..];
                Some(read_slice(&mut src))
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
