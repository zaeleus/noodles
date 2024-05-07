mod alternative_allele;
mod contig;
mod filter;
mod format;
mod info;
mod meta;
mod other;

use std::io::{self, Write};

pub(crate) use self::{
    alternative_allele::write_alternative_allele, contig::write_contig, filter::write_filter,
    format::write_format, info::write_info, meta::write_meta, other::write_other,
};
use crate::{
    header::record::value::map::{self, tag, OtherFields},
    io::writer::header::record::write_separator,
};

pub(crate) fn write_map<W, I, F>(writer: &mut W, id: I, f: F) -> io::Result<()>
where
    W: Write,
    I: AsRef<str>,
    F: Fn(&mut W) -> io::Result<()>,
{
    write_prefix(writer)?;
    write_id_field(writer, id)?;
    f(writer)?;
    write_suffix(writer)?;
    Ok(())
}

pub(crate) fn write_other_map<W, I, F>(
    writer: &mut W,
    id_tag: &map::other::Tag,
    id: I,
    f: F,
) -> io::Result<()>
where
    W: Write,
    I: AsRef<str>,
    F: Fn(&mut W) -> io::Result<()>,
{
    write_prefix(writer)?;
    write_value_field(writer, id_tag, id)?;
    f(writer)?;
    write_suffix(writer)?;
    Ok(())
}

fn write_id_field<W, I>(writer: &mut W, id: I) -> io::Result<()>
where
    W: Write,
    I: AsRef<str>,
{
    use crate::header::record::value::map::tag::ID;

    write_value_field(writer, ID, id)?;

    Ok(())
}

fn write_type_field<W, T>(writer: &mut W, ty: T) -> io::Result<()>
where
    W: Write,
    T: AsRef<str>,
{
    use crate::header::record::value::map::tag::TYPE;

    write_delimiter(writer)?;
    write_value_field(writer, TYPE, ty)?;

    Ok(())
}

fn write_description_field<W>(writer: &mut W, description: &str) -> io::Result<()>
where
    W: Write,
{
    use crate::header::record::value::map::tag::DESCRIPTION;

    write_delimiter(writer)?;
    write_string_field(writer, DESCRIPTION, description)?;

    Ok(())
}

fn write_other_fields<W, S>(writer: &mut W, other_fields: &OtherFields<S>) -> io::Result<()>
where
    W: Write,
    S: tag::Standard,
{
    for (key, value) in other_fields {
        write_delimiter(writer)?;
        write_string_field(writer, key, value)?;
    }

    Ok(())
}

fn write_value_field<W, K, V>(writer: &mut W, key: K, value: V) -> io::Result<()>
where
    W: Write,
    K: AsRef<str>,
    V: AsRef<str>,
{
    write_key(writer, key)?;
    write_separator(writer)?;
    write_value(writer, value)?;
    Ok(())
}

fn write_string_field<W, K, V>(writer: &mut W, key: K, value: V) -> io::Result<()>
where
    W: Write,
    K: AsRef<str>,
    V: AsRef<str>,
{
    write_key(writer, key)?;
    write_separator(writer)?;
    write_escaped_string(writer, value)?;
    Ok(())
}

fn write_key<W, K>(writer: &mut W, key: K) -> io::Result<()>
where
    W: Write,
    K: AsRef<str>,
{
    writer.write_all(key.as_ref().as_bytes())
}

fn write_value<W, V>(writer: &mut W, value: V) -> io::Result<()>
where
    W: Write,
    V: AsRef<str>,
{
    writer.write_all(value.as_ref().as_bytes())
}

fn write_prefix<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const PREFIX: u8 = b'<';
    writer.write_all(&[PREFIX])
}

fn write_suffix<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const SUFFIX: u8 = b'>';
    writer.write_all(&[SUFFIX])
}

fn write_delimiter<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: u8 = b',';
    writer.write_all(&[DELIMITER])
}

fn write_escaped_string<W, V>(writer: &mut W, s: V) -> io::Result<()>
where
    W: Write,
    V: AsRef<str>,
{
    write!(writer, "\"")?;

    for c in s.as_ref().chars() {
        if matches!(c, '"' | '\\') {
            write!(writer, "\\")?;
        }

        write!(writer, "{c}")?;
    }

    write!(writer, "\"")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_escaped_string() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, s: &str, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_escaped_string(buf, s)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, r#"noodles"#, br#""noodles""#)?;
        t(&mut buf, "noodles=🍜", b"\"noodles=\xf0\x9f\x8d\x9c\"")?;
        t(&mut buf, r#"noodles-"vcf""#, br#""noodles-\"vcf\"""#)?;
        t(&mut buf, r"noodles\vcf", br#""noodles\\vcf""#)?;

        Ok(())
    }
}
