mod byte;
mod byte_array;
mod integer;

use std::io;

pub use self::{byte::Byte, byte_array::ByteArray, integer::Integer};
use crate::io::{
    reader::container::slice::records::ExternalDataReaders,
    writer::container::slice::records::ExternalDataWriters, BitReader, BitWriter,
};

pub trait Decode<'de> {
    type Value;

    fn decode(
        &self,
        core_data_reader: &mut BitReader<'de>,
        external_data_readers: &mut ExternalDataReaders<'de>,
    ) -> io::Result<Self::Value>;
}

pub trait Encode<'en> {
    type Value;

    fn encode(
        &self,
        core_data_writer: &mut BitWriter,
        external_data_writers: &mut ExternalDataWriters,
        value: Self::Value,
    ) -> io::Result<()>;
}
