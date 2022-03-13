mod field;

use std::io;

use bytes::BufMut;
use noodles_sam::record::Data;

pub use self::field::put_field;

pub fn put_data<B>(dst: &mut B, data: &Data) -> io::Result<()>
where
    B: BufMut,
{
    for field in data.values() {
        put_field(dst, field)?;
    }

    Ok(())
}
