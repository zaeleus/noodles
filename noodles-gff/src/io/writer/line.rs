mod comment;
mod directive;
mod record;

use std::io::{self, Write};

use self::comment::write_comment;
pub(super) use self::{directive::write_directive, record::write_record};
use crate::LineBuf;

pub(super) fn write_line<W>(writer: &mut W, line: &LineBuf) -> io::Result<()>
where
    W: Write,
{
    match line {
        LineBuf::Directive(directive) => write_directive(writer, directive)?,
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

#[cfg(test)]
mod tests {
    use bstr::BString;

    use super::*;
    use crate::{
        directive_buf::{key, Value},
        feature::RecordBuf,
        DirectiveBuf,
    };

    #[test]
    fn test_write_line() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, line: &LineBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_line(buf, line)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let line = LineBuf::Directive(DirectiveBuf::new(
            key::GFF_VERSION,
            Some(Value::GffVersion(Default::default())),
        ));
        t(&mut buf, &line, b"##gff-version 3\n")?;

        let line = LineBuf::Comment(BString::from("noodles"));
        t(&mut buf, &line, b"#noodles\n")?;

        let line = LineBuf::Record(RecordBuf::default());
        t(&mut buf, &line, b".\t.\t.\t1\t1\t.\t.\t.\t.\n")?;

        Ok(())
    }
}
