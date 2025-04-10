mod comment;
mod record;

use std::io::{self, Write};

use self::comment::write_comment;
pub(super) use self::record::write_record;
use crate::LineBuf;

pub(super) fn write_line<W>(writer: &mut W, line: &LineBuf) -> io::Result<()>
where
    W: Write,
{
    match line {
        LineBuf::Comment(s) => write_comment(writer, s.as_ref())?,
        LineBuf::Record(record) => write_record(writer, record)?,
    }

    write_newline(writer)?;

    Ok(())
}

pub(super) fn write_newline<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const LINE_FEED: u8 = b'\n';
    writer.write_all(&[LINE_FEED])
}
