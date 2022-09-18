#![allow(dead_code)]

mod indexed_raw_reader;
mod raw_reader;

pub use self::{indexed_raw_reader::IndexedRawReader, raw_reader::RawReader};

pub enum Inner<R> {
    Raw(RawReader<R>),
    IndexedRaw(IndexedRawReader<R>),
}
