mod compression_method;
mod content_type;

pub use self::{compression_method::CompressionMethod, content_type::ContentType};

use std::{borrow::Cow, io::Read};

use bzip2::read::BzDecoder;
use flate2::read::GzDecoder;
use xz2::read::XzDecoder;

use crate::{num::Itf8, rans::rans_decode};

#[derive(Debug, Default)]
pub struct Block {
    compression_method: CompressionMethod,
    content_type: u8,
    content_id: Itf8,
    uncompressed_len: Itf8,
    data: Vec<u8>,
    crc32: [u8; 4],
}

impl Block {
    pub fn compression_method(&self) -> CompressionMethod {
        self.compression_method
    }

    pub fn compression_method_mut(&mut self) -> &mut CompressionMethod {
        &mut self.compression_method
    }

    pub fn content_type(&self) -> u8 {
        self.content_type
    }

    pub fn content_type_mut(&mut self) -> &mut u8 {
        &mut self.content_type
    }

    pub fn content_id(&self) -> Itf8 {
        self.content_id
    }

    pub fn content_id_mut(&mut self) -> &mut Itf8 {
        &mut self.content_id
    }

    pub fn uncompressed_len(&self) -> Itf8 {
        self.uncompressed_len
    }

    pub fn uncompressed_len_mut(&mut self) -> &mut Itf8 {
        &mut self.uncompressed_len
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut Vec<u8> {
        &mut self.data
    }

    pub fn decompressed_data(&self) -> Cow<[u8]> {
        match self.compression_method {
            CompressionMethod::None => Cow::from(self.data()),
            CompressionMethod::Gzip => {
                let mut reader = GzDecoder::new(self.data());
                let mut buf = Vec::with_capacity(self.uncompressed_len as usize);
                reader.read_to_end(&mut buf).expect("invalid gzip data");
                Cow::from(buf)
            }
            CompressionMethod::Bzip2 => {
                let mut reader = BzDecoder::new(self.data());
                let mut buf = Vec::with_capacity(self.uncompressed_len as usize);
                reader.read_to_end(&mut buf).expect("invalid bzip2 data");
                Cow::from(buf)
            }
            CompressionMethod::Lzma => {
                let mut reader = XzDecoder::new(self.data());
                let mut buf = Vec::with_capacity(self.uncompressed_len as usize);
                reader.read_to_end(&mut buf).expect("invalid lzma data");
                Cow::from(buf)
            }
            CompressionMethod::Rans => {
                let mut buf = self.data();
                rans_decode(&mut buf)
                    .map(Cow::from)
                    .expect("invalid rans data")
            }
        }
    }

    pub fn crc32(&self) -> &[u8; 4] {
        &self.crc32
    }

    pub fn crc32_mut(&mut self) -> &mut [u8; 4] {
        &mut self.crc32
    }
}
