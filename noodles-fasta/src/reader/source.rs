use std::{
    fs::File,
    io::{self, BufRead, BufReader, Read, Seek, SeekFrom},
};

use noodles_bgzf::{self as bgzf, gzi};

pub enum Source<R> {
    File(BufReader<R>),
    Bgzip(bgzf::Reader<R>),
    IndexedBgzip(bgzf::Reader<R>, gzi::Index),
}

impl<R> Read for Source<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::File(inner) => inner.read(buf),
            Self::Bgzip(inner) => inner.read(buf),
            Self::IndexedBgzip(inner, _) => inner.read(buf),
        }
    }
}

impl<R> BufRead for Source<R>
where
    R: Read,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            Self::File(inner) => inner.fill_buf(),
            Self::Bgzip(inner) => inner.fill_buf(),
            Self::IndexedBgzip(inner, _) => inner.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            Self::File(inner) => inner.consume(amt),
            Self::Bgzip(inner) => inner.consume(amt),
            Self::IndexedBgzip(inner, _) => inner.consume(amt),
        }
    }
}

impl<R> Seek for Source<R>
where
    R: Read + Seek,
{
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        match self {
            Self::File(inner) => inner.seek(pos),
            Self::Bgzip(_) => Err(io::Error::from(io::ErrorKind::Unsupported)),
            Self::IndexedBgzip(inner, index) => match pos {
                SeekFrom::Start(p) => inner.seek_by_uncompressed_position(index, p),
                _ => Err(io::Error::from(io::ErrorKind::Unsupported)),
            },
        }
    }
}

impl From<BufReader<File>> for Source<File> {
    fn from(reader: BufReader<File>) -> Self {
        Self::File(reader)
    }
}

impl<R> From<bgzf::Reader<R>> for Source<R>
where
    R: Read,
{
    fn from(reader: bgzf::Reader<R>) -> Self {
        Self::Bgzip(reader)
    }
}
