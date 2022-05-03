//! BAM record data field reader.

pub mod field;

use bytes::BytesMut;

pub(crate) use self::field::get_field;
use crate::record::Data;

pub(crate) fn get_data(src: &mut BytesMut, data: &mut Data) {
    data.clear();
    data.buf = src.split();
}
