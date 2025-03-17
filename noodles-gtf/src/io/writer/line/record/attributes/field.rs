use std::io::{self, Write};

use noodles_gff::feature::record_buf::attributes::field::Value;

use super::write_separator;

pub(super) fn write_field<W>(writer: &mut W, key: &str, value: &Value) -> io::Result<()>
where
    W: Write,
{
    for (i, v) in value.iter().enumerate() {
        if i > 0 {
            write_separator(writer)?;
        }

        writer.write_all(key.as_bytes())?;
        write_separator(writer)?;
        write!(writer, r#""{}""#, v)?;
        write_terminator(writer)?;
    }

    Ok(())
}

fn write_terminator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const TERMINATOR: u8 = b';';
    writer.write_all(&[TERMINATOR])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_field() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, key: &str, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_field(buf, key, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, "gene_id", &Value::from("g0"), br#"gene_id "g0";"#)?;

        t(
            &mut buf,
            "tag",
            &Value::from(vec![String::from("nd"), String::from("ls")]),
            br#"tag "nd"; tag "ls";"#,
        )?;

        Ok(())
    }
}
