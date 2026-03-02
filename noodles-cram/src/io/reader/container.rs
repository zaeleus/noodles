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
use crate::{
    container::{CompressionHeader, Header},
    file_definition::Version,
};

/// A CRAM container.
#[derive(Default)]
pub struct Container {
    pub(crate) header: Header,
    pub(crate) src: Vec<u8>,
    pub(crate) version: Version,
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

        read_compression_header(&mut src, self.version)
    }

    /// Returns the iterator over slices.
    pub fn slices(&self) -> impl Iterator<Item = io::Result<Slice<'_>>> + '_ {
        let landmarks = &self.header.landmarks;
        let version = self.version;
        let mut i = 0;

        iter::from_fn(move || {
            if landmarks.is_empty() {
                return None;
            }

            if i < landmarks.len() - 1 {
                let (start, end) = (landmarks[i], landmarks[i + 1]);
                i += 1;
                let mut src = &self.src[start..end];
                Some(read_slice(&mut src, version))
            } else if i < landmarks.len() {
                let start = landmarks[i];
                i += 1;
                let mut src = &self.src[start..];
                Some(read_slice(&mut src, version))
            } else {
                None
            }
        })
    }

    /// Returns slices decoded in parallel.
    ///
    /// This is equivalent to [`Self::slices`] but uses rayon to decode slices concurrently.
    #[cfg(feature = "parallel")]
    pub fn slices_par(&self) -> io::Result<Vec<Slice<'_>>> {
        use rayon::prelude::*;

        let landmarks = &self.header.landmarks;
        let version = self.version;
        let src = &self.src;

        (0..landmarks.len())
            .into_par_iter()
            .map(|i| {
                let start = landmarks[i];
                let end = landmarks.get(i + 1).copied().unwrap_or(src.len());
                let mut slice_src = &src[start..end];
                read_slice(&mut slice_src, version)
            })
            .collect()
    }
}

pub fn read_container<R>(
    reader: &mut R,
    container: &mut Container,
    version: Version,
) -> io::Result<usize>
where
    R: Read,
{
    container.version = version;

    match read_header(reader, &mut container.header, version)? {
        0 => Ok(0),
        len => {
            container.src.resize(len, 0);
            reader.read_exact(&mut container.src)?;
            Ok(len)
        }
    }
}
