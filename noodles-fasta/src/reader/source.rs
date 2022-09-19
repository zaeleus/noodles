use std::{
    fs::File,
    io::{self, Read},
};

use noodles_bgzf as bgzf;

pub enum Source<R> {
    File(R),
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

impl From<File> for Source<File> {
    fn from(file: File) -> Self {
        Self::File(file)
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
