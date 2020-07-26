mod compression_method;
mod content_type;

pub use self::{compression_method::CompressionMethod, content_type::ContentType};

use std::{
    borrow::Cow,
    io::{self, Read},
    mem,
};

use bzip2::read::BzDecoder;
use flate2::read::GzDecoder;
use xz2::read::XzDecoder;

use crate::{
    num::{write_itf8, Itf8},
    rans::rans_decode,
};

#[derive(Clone, Debug)]
pub struct Block {
    compression_method: CompressionMethod,
    content_type: u8,
    content_id: Itf8,
    uncompressed_len: Itf8,
    data: Vec<u8>,
    crc32: u32,
}

impl Block {
    pub fn new(
        compression_method: CompressionMethod,
        content_type: u8,
        content_id: Itf8,
        uncompressed_len: Itf8,
        data: Vec<u8>,
        crc32: u32,
    ) -> Self {
        Self {
            compression_method,
            content_type,
            content_id,
            uncompressed_len,
            data,
            crc32,
        }
    }

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

    pub fn crc32(&self) -> u32 {
        self.crc32
    }

    pub fn crc32_mut(&mut self) -> &mut u32 {
        &mut self.crc32
    }

    pub fn len(&self) -> io::Result<usize> {
        // method (1) + block content type ID (1)
        let mut len = 2 * mem::size_of::<u8>();

        let mut buf = Vec::new();
        write_itf8(&mut buf, self.content_id())?;
        len += buf.len();

        buf.clear();
        write_itf8(&mut buf, self.data().len() as i32)?;
        len += buf.len();

        buf.clear();
        write_itf8(&mut buf, self.uncompressed_len())?;
        len += buf.len();

        len += self.data.len();

        // crc32 (4)
        len += mem::size_of::<u32>();

        Ok(len)
    }
}
