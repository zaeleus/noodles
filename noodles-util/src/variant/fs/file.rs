#![expect(dead_code)]

use noodles_bgzf as bgzf;

use crate::variant::io::{CompressionMethod, reader::builder::detect_compression_method};

use std::{
    io::{self, BufRead, BufReader, Read},
    path::Path,
};

pub enum File {
    Bgzf(bgzf::io::Reader<BufReader<std::fs::File>>),
    Raw(BufReader<std::fs::File>),
}

impl File {
    pub fn open<P>(src: P) -> io::Result<Self>
    where
        P: AsRef<Path>,
    {
        let mut reader = std::fs::File::open(src).map(BufReader::new)?;

        match detect_compression_method(&mut reader)? {
            None => Ok(Self::Raw(reader)),
            Some(CompressionMethod::Bgzf) => Ok(Self::Bgzf(bgzf::io::Reader::new(reader))),
        }
    }
}

impl Read for File {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::Bgzf(reader) => reader.read(buf),
            Self::Raw(reader) => reader.read(buf),
        }
    }
}

impl BufRead for File {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            Self::Bgzf(reader) => reader.fill_buf(),
            Self::Raw(reader) => reader.fill_buf(),
        }
    }

    fn consume(&mut self, amount: usize) {
        match self {
            Self::Bgzf(reader) => reader.consume(amount),
            Self::Raw(reader) => reader.consume(amount),
        }
    }
}
