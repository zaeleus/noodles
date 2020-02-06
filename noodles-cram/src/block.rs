use std::{convert::TryFrom, io::Read};

use flate2::read::GzDecoder;

use crate::num::Itf8;

#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub enum CompressionMethod {
    None,
    Gzip,
    Bzip2,
    Lzma,
    Rans,
}

impl TryFrom<u8> for CompressionMethod {
    type Error = ();

    fn try_from(b: u8) -> Result<Self, Self::Error> {
        match b {
            0 => Ok(Self::None),
            1 => Ok(Self::Gzip),
            2 => Ok(Self::Bzip2),
            3 => Ok(Self::Lzma),
            4 => Ok(Self::Rans),
            _ => Err(()),
        }
    }
}

#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub enum ContentType {
    FileHeader,
    CompressionHeader,
    SliceHeader,
    Reserved,
    ExternalData,
    CoreData,
}

impl TryFrom<u8> for ContentType {
    type Error = ();

    fn try_from(b: u8) -> Result<Self, Self::Error> {
        match b {
            0 => Ok(Self::FileHeader),
            1 => Ok(Self::CompressionHeader),
            2 => Ok(Self::SliceHeader),
            3 => Ok(Self::Reserved),
            4 => Ok(Self::ExternalData),
            5 => Ok(Self::CoreData),
            _ => Err(()),
        }
    }
}

#[derive(Debug, Default)]
pub struct Block {
    compression_method: u8,
    content_type: u8,
    content_id: Itf8,
    data: Vec<u8>,
    crc32: [u8; 4],
}

impl Block {
    pub fn compression_method(&self) -> u8 {
        self.compression_method
    }

    pub fn compression_method_mut(&mut self) -> &mut u8 {
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

    pub fn data(&self) -> &[u8] {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut Vec<u8> {
        &mut self.data
    }

    pub fn decompressed_data(&self) -> Vec<u8> {
        let compression_method = CompressionMethod::try_from(self.compression_method())
            .expect("invalid compression method");

        match compression_method {
            CompressionMethod::None => self.data().to_vec(),
            CompressionMethod::Gzip => {
                let mut reader = GzDecoder::new(self.data());
                let mut buf = Vec::new();
                reader.read_to_end(&mut buf).expect("invalid gzip data");
                buf
            }
            _ => todo!(),
        }
    }
}
