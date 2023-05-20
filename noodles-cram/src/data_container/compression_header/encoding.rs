pub mod codec;
mod kind;

use bytes::Buf;

pub use self::kind::Kind;

use std::io;

use crate::{io::BitReader, reader::record::ExternalDataReaders};

pub trait Decode {
    type Value;

    fn decode<R, S>(
        &self,
        core_data_reader: &mut BitReader<R>,
        external_data_readers: &mut ExternalDataReaders<S>,
    ) -> io::Result<Self::Value>
    where
        R: Buf,
        S: Buf;
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Encoding<C>(C);

impl<C> Encoding<C> {
    pub fn new(codec: C) -> Self {
        Self(codec)
    }

    pub fn get(&self) -> &C {
        &self.0
    }
}

impl<C> Encoding<C>
where
    C: Decode,
{
    pub fn decode<R, S>(
        &self,
        core_data_reader: &mut BitReader<R>,
        external_data_readers: &mut ExternalDataReaders<S>,
    ) -> io::Result<C::Value>
    where
        R: Buf,
        S: Buf,
    {
        self.get().decode(core_data_reader, external_data_readers)
    }
}
