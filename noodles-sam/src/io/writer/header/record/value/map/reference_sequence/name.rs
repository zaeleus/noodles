use std::io::{self, Write};

use crate::header::record::value::map::reference_sequence::tag;

pub(super) fn write_name_field<W>(writer: &mut W, name: &[u8]) -> io::Result<()>
where
    W: Write,
{
    use crate::io::writer::header::record::{value::map::write_separator, write_delimiter};

    write_delimiter(writer)?;
    writer.write_all(tag::NAME.as_ref())?;
    write_separator(writer)?;
    write_value(writer, name)?;

    Ok(())
}

fn write_value<W>(writer: &mut W, name: &[u8]) -> io::Result<()>
where
    W: Write,
{
    if is_valid_name(name) {
        writer.write_all(name)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid reference sequence name",
        ))
    }
}

//  § 1.2.1 "Character set restrictions" (2023-05-24): "...`[:rname:∧*=][:rname:]*`."
fn is_valid_name(name: &[u8]) -> bool {
    let mut iter = name.iter().copied();

    if let Some(b) = iter.next() {
        if b == b'*' || b == b'=' || !is_valid_name_char(b) {
            return false;
        }

        iter.all(is_valid_name_char)
    } else {
        false
    }
}

fn is_valid_name_char(b: u8) -> bool {
    b.is_ascii_graphic()
        && !matches!(
            b,
            b'\\'
                | b','
                | b'"'
                | b'`'
                | b'\''
                | b'('
                | b')'
                | b'['
                | b']'
                | b'{'
                | b'}'
                | b'<'
                | b'>',
        )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_valid_name() {
        assert!(is_valid_name(b"sq0"));
        assert!(is_valid_name(b"sq0*"));
        assert!(is_valid_name(b"sq0="));

        assert!(!is_valid_name(b""));
        assert!(!is_valid_name(b"sq 0"));
        assert!(!is_valid_name(b"sq[0]"));
        assert!(!is_valid_name(b">sq0"));
        assert!(!is_valid_name(b"*sq0"));
        assert!(!is_valid_name(b"=sq0"));
    }
}
