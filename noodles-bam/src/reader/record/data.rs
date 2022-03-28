//! BAM record data field reader.

pub mod field;

pub(crate) use self::field::get_field;

use std::io;

use bytes::Buf;
use noodles_sam as sam;

pub(super) fn get_data<B>(src: &mut B, data: &mut sam::record::Data) -> io::Result<()>
where
    B: Buf,
{
    data.clear();

    while let Some(field) = get_field(src)? {
        data.insert(field);
    }

    Ok(())
}
