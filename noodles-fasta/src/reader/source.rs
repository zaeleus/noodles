use std::{
    fs::File,
    io::{self, BufRead, BufReader, Read},
};

use noodles_bgzf as bgzf;

pub enum Source<R> {
    File(BufReader<R>),
    Bgzip(bgzf::Reader<R>),
}

impl<R> Read for Source<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::File(inner) => inner.read(buf),
            Self::Bgzip(inner) => inner.read(buf),
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
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            Self::File(inner) => inner.consume(amt),
            Self::Bgzip(inner) => inner.consume(amt),
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
