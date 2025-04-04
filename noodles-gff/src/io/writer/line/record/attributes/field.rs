mod tag;
mod value;

use std::{
    borrow::Cow,
    io::{self, Write},
};

use bstr::BStr;
use percent_encoding::{AsciiSet, CONTROLS};

use self::{tag::write_tag, value::write_value};
use crate::feature::record::attributes::field::Value;

pub(super) fn write_field<W>(writer: &mut W, key: &BStr, value: &Value) -> io::Result<()>
where
    W: Write,
{
    write_tag(writer, key)?;
    write_separator(writer)?;
    write_value(writer, value)?;
    Ok(())
}

fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: u8 = b'=';
    writer.write_all(&[SEPARATOR])
}

fn percent_encode(s: &BStr) -> Cow<'_, str> {
    const PERCENT_ENCODE_SET: &AsciiSet = &CONTROLS
        .add(b'\t')
        .add(b'\n')
        .add(b'\r')
        .add(b'%')
        .add(b';')
        .add(b'=')
        .add(b'&')
        .add(b',');

    percent_encoding::percent_encode(s, PERCENT_ENCODE_SET).into()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_field() -> io::Result<()> {
        let mut buf = Vec::new();

        write_field(
            &mut buf,
            BStr::new("ID"),
            &Value::String(Cow::from(BStr::new("0"))),
        )?;

        assert_eq!(buf, b"ID=0");

        Ok(())
    }
}
