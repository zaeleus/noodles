use std::{
    borrow::Cow,
    io::{self, Write},
};

use percent_encoding::{utf8_percent_encode, AsciiSet, NON_ALPHANUMERIC};

pub(super) fn write_reference_sequence_name<W>(
    writer: &mut W,
    reference_sequence_name: &str,
) -> io::Result<()>
where
    W: Write,
{
    let s = percent_encode(reference_sequence_name);
    writer.write_all(s.as_bytes())
}

// § "Column 1: 'seqid'" (2020-08-18): "IDs may contain any characters, but must escape any
// characters not in the set [a-zA-Z0-9.:^*$@!+_?-|]."
fn percent_encode(s: &str) -> Cow<'_, str> {
    const PERCENT_ENCODE_SET: &AsciiSet = &NON_ALPHANUMERIC
        .remove(b'.')
        .remove(b':')
        .remove(b'^')
        .remove(b'*')
        .remove(b'$')
        .remove(b'@')
        .remove(b'!')
        .remove(b'+')
        .remove(b'_')
        .remove(b'?')
        .remove(b'-')
        .remove(b'|');

    utf8_percent_encode(s, PERCENT_ENCODE_SET).into()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_reference_sequence_name() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, reference_sequence_name: &str, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_reference_sequence_name(buf, reference_sequence_name)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, "sq0", b"sq0")?;
        t(&mut buf, "sq0.:^*$@!+_?-|", b"sq0.:^*$@!+_?-|")?;
        t(&mut buf, "sq 0", b"sq%200")?;
        t(&mut buf, ">sq0", b"%3Esq0")?;

        Ok(())
    }
}
