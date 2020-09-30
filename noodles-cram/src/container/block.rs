mod builder;
mod compression_method;
mod content_type;

pub use self::{
    builder::Builder, compression_method::CompressionMethod, content_type::ContentType,
};

use std::{
    borrow::Cow,
    io::{self, Read},
    mem,
};

use bzip2::read::BzDecoder;
use flate2::read::GzDecoder;
use xz2::read::XzDecoder;

use crate::{
    num::{itf8, Itf8},
    rans::rans_decode,
};

// § 9 End of file container (2020-06-22)
const EOF_DATA: [u8; 6] = [0x01, 0x00, 0x01, 0x00, 0x01, 0x00];
const EOF_CRC32: u32 = 0x4b_01_63_ee;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Block {
    compression_method: CompressionMethod,
    content_type: ContentType,
    content_id: Itf8,
    uncompressed_len: Itf8,
    data: Vec<u8>,
    crc32: u32,
}

#[allow(clippy::len_without_is_empty)]
impl Block {
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Creates a block used in the EOF container.
    pub fn eof() -> Self {
        Self::builder()
            .set_content_type(ContentType::CompressionHeader)
            .set_uncompressed_len(EOF_DATA.len() as Itf8)
            .set_data(EOF_DATA.to_vec())
            .set_crc32(EOF_CRC32)
            .build()
    }

    pub fn compression_method(&self) -> CompressionMethod {
        self.compression_method
    }

    pub fn content_type(&self) -> ContentType {
        self.content_type
    }

    pub fn content_id(&self) -> Itf8 {
        self.content_id
    }

    pub fn uncompressed_len(&self) -> Itf8 {
        self.uncompressed_len
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }

    pub fn decompressed_data(&self) -> io::Result<Cow<'_, [u8]>> {
        match self.compression_method {
            CompressionMethod::None => Ok(Cow::from(self.data())),
            CompressionMethod::Gzip => {
                let mut reader = GzDecoder::new(self.data());
                let mut buf = Vec::with_capacity(self.uncompressed_len as usize);
                reader.read_to_end(&mut buf)?;
                Ok(Cow::from(buf))
            }
            CompressionMethod::Bzip2 => {
                let mut reader = BzDecoder::new(self.data());
                let mut buf = Vec::with_capacity(self.uncompressed_len as usize);
                reader.read_to_end(&mut buf)?;
                Ok(Cow::from(buf))
            }
            CompressionMethod::Lzma => {
                let mut reader = XzDecoder::new(self.data());
                let mut buf = Vec::with_capacity(self.uncompressed_len as usize);
                reader.read_to_end(&mut buf)?;
                Ok(Cow::from(buf))
            }
            CompressionMethod::Rans => {
                let mut buf = self.data();
                rans_decode(&mut buf).map(Cow::from)
            }
        }
    }

    pub fn crc32(&self) -> u32 {
        self.crc32
    }

    pub fn len(&self) -> usize {
        // method
        mem::size_of::<u8>()
            // block content type ID
            + mem::size_of::<u8>()
            + itf8::size_of(self.content_id())
            + itf8::size_of(self.data.len() as Itf8)
            + itf8::size_of(self.uncompressed_len())
            + self.data.len()
            // crc32
            + mem::size_of::<u32>()
    }
}
