mod definition;

pub(crate) use self::definition::read_definition;

use std::io::{self, BufRead, Read};

use crate::Record;

const LINE_FEED: u8 = b'\n';
const CARRIAGE_RETURN: u8 = b'\r';

pub(super) fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: BufRead,
{
    record.clear();

    // Fast path: try to read the entire record from the current buffer.
    if let Some(n) = try_read_record_bulk(reader, record)? {
        return Ok(n);
    }

    // Slow path: record spans the buffer boundary. Fall back to line-by-line.
    // The bulk path only called fill_buf() (not consume()), so the cursor
    // has not advanced. The bulk path returns None before writing any fields,
    // so the record is still in the cleared state from above.
    let mut len = match read_definition(reader, record.definition_mut()) {
        Ok(0) => return Ok(0),
        Ok(n) => n,
        Err(e) => return Err(e),
    };

    len += read_line(reader, record.sequence_mut())?;
    len += consume_plus_line(reader)?;
    len += read_line(reader, record.quality_scores_mut())?;

    Ok(len)
}

/// Attempt to read a complete FASTQ record from the current BufRead buffer
/// in a single pass. Returns `None` if the record spans the buffer boundary,
/// in which case the caller should fall back to line-by-line parsing.
///
/// # Preconditions
///
/// The BufRead cursor must be positioned at the start of a record (i.e.,
/// the next byte is `@` or EOF). This is guaranteed because `read_record`
/// always consumes a complete record before calling this again.
fn try_read_record_bulk<R>(reader: &mut R, record: &mut Record) -> io::Result<Option<usize>>
where
    R: BufRead,
{
    use memchr::memchr;

    let src = reader.fill_buf()?;

    if src.is_empty() {
        return Ok(Some(0));
    }

    if src[0] != b'@' {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid name prefix",
        ));
    }

    // Find 4 newlines in the buffer.
    let mut newlines = [0usize; 4];
    let mut pos = 0;
    for nl in &mut newlines {
        match memchr(LINE_FEED, &src[pos..]) {
            Some(i) => {
                *nl = pos + i;
                pos = pos + i + 1;
            }
            None => return Ok(None), // record spans buffer boundary
        }
    }

    // Validate '+' prefix on line 3.
    let plus_start = newlines[1] + 1;
    if plus_start >= src.len() || src[plus_start] != b'+' {
        return Ok(None); // malformed or spans boundary; let fallback produce the error
    }

    let total_len = newlines[3] + 1;

    // Line 1: @name [description]
    let def_line = &src[1..newlines[0]];
    let def_line = def_line.strip_suffix(b"\r").unwrap_or(def_line);
    if let Some(i) = memchr::memchr2(b' ', b'\t', def_line) {
        record.definition_mut().name_mut().extend(&def_line[..i]);
        record
            .definition_mut()
            .description_mut()
            .extend(&def_line[i + 1..]);
    } else {
        record.definition_mut().name_mut().extend(def_line);
    }

    // Line 2: sequence
    let seq = &src[newlines[0] + 1..newlines[1]];
    let seq = seq.strip_suffix(b"\r").unwrap_or(seq);
    record.sequence_mut().extend_from_slice(seq);

    // Line 3: + (skipped)
    // Line 4: quality
    let qual = &src[newlines[2] + 1..newlines[3]];
    let qual = qual.strip_suffix(b"\r").unwrap_or(qual);
    record.quality_scores_mut().extend_from_slice(qual);

    reader.consume(total_len);

    Ok(Some(total_len))
}

fn read_line<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    match reader.read_until(LINE_FEED, buf) {
        Ok(0) => Ok(0),
        Ok(n) => {
            if buf.ends_with(&[LINE_FEED]) {
                buf.pop();

                if buf.ends_with(&[CARRIAGE_RETURN]) {
                    buf.pop();
                }
            }

            Ok(n)
        }
        Err(e) => Err(e),
    }
}

fn consume_line<R>(reader: &mut R) -> io::Result<usize>
where
    R: BufRead,
{
    use memchr::memchr;

    let mut is_eol = false;
    let mut len = 0;

    loop {
        let src = reader.fill_buf()?;

        if src.is_empty() || is_eol {
            break;
        }

        let n = match memchr(LINE_FEED, src) {
            Some(i) => {
                is_eol = true;
                i + 1
            }
            None => src.len(),
        };

        reader.consume(n);

        len += n;
    }

    Ok(len)
}

fn read_u8<R>(reader: &mut R) -> io::Result<u8>
where
    R: Read,
{
    let mut buf = [0; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn consume_plus_line<R>(reader: &mut R) -> io::Result<usize>
where
    R: BufRead,
{
    const PREFIX: u8 = b'+';

    match read_u8(reader)? {
        PREFIX => consume_line(reader).map(|n| n + 1),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid description prefix",
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_line() -> io::Result<()> {
        let mut buf = Vec::new();

        let data = b"noodles\n";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, b"noodles");

        let data = b"noodles\r\n";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, b"noodles");

        let data = b"noodles";
        let mut reader = &data[..];
        buf.clear();
        read_line(&mut reader, &mut buf)?;
        assert_eq!(buf, b"noodles");

        Ok(())
    }

    #[test]
    fn test_consume_plus_line() -> io::Result<()> {
        let data = b"+r0\n";
        let mut reader = &data[..];
        consume_plus_line(&mut reader)?;

        let data = b"r0\n";
        let mut reader = &data[..];
        assert!(matches!(
            consume_plus_line(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_read_record_bulk() -> io::Result<()> {
        let data = b"@r0 desc0\nACGT\n+\nNDLS\n@r1\nTGCA\n+\nSLDN\n";
        let mut reader = &data[..];
        let mut record = Record::default();

        read_record(&mut reader, &mut record)?;
        assert_eq!(record.name(), &b"r0"[..]);
        assert_eq!(record.description(), &b"desc0"[..]);
        assert_eq!(record.sequence(), b"ACGT");
        assert_eq!(record.quality_scores(), b"NDLS");

        read_record(&mut reader, &mut record)?;
        assert_eq!(record.name(), &b"r1"[..]);
        assert_eq!(record.description(), &b""[..]);
        assert_eq!(record.sequence(), b"TGCA");
        assert_eq!(record.quality_scores(), b"SLDN");

        assert_eq!(read_record(&mut reader, &mut record)?, 0);
        Ok(())
    }

    #[test]
    fn test_read_record_crlf() -> io::Result<()> {
        let data = b"@r0\r\nACGT\r\n+\r\nNDLS\r\n";
        let mut reader = &data[..];
        let mut record = Record::default();

        read_record(&mut reader, &mut record)?;
        assert_eq!(record.name(), &b"r0"[..]);
        assert_eq!(record.sequence(), b"ACGT");
        assert_eq!(record.quality_scores(), b"NDLS");
        Ok(())
    }

    #[test]
    fn test_read_record_with_plus_description() -> io::Result<()> {
        let data = b"@r0\nACGT\n+r0\nNDLS\n";
        let mut reader = &data[..];
        let mut record = Record::default();

        read_record(&mut reader, &mut record)?;
        assert_eq!(record.name(), &b"r0"[..]);
        assert_eq!(record.sequence(), b"ACGT");
        assert_eq!(record.quality_scores(), b"NDLS");
        Ok(())
    }

    #[test]
    fn test_read_record_tab_description() -> io::Result<()> {
        let data = b"@r0\tLN:4\nACGT\n+\nNDLS\n";
        let mut reader = &data[..];
        let mut record = Record::default();

        read_record(&mut reader, &mut record)?;
        assert_eq!(record.name(), &b"r0"[..]);
        assert_eq!(record.description(), &b"LN:4"[..]);
        assert_eq!(record.sequence(), b"ACGT");
        assert_eq!(record.quality_scores(), b"NDLS");
        Ok(())
    }

    #[test]
    fn test_read_record_empty_sequence() -> io::Result<()> {
        let data = b"@r0\n\n+\n\n";
        let mut reader = &data[..];
        let mut record = Record::default();

        read_record(&mut reader, &mut record)?;
        assert_eq!(record.name(), &b"r0"[..]);
        assert_eq!(record.sequence(), b"");
        assert_eq!(record.quality_scores(), b"");
        Ok(())
    }

    #[test]
    fn test_read_record_quality_contains_at() -> io::Result<()> {
        let data = b"@r0\nACGT\n+\n@DLS\n";
        let mut reader = &data[..];
        let mut record = Record::default();

        read_record(&mut reader, &mut record)?;
        assert_eq!(record.sequence(), b"ACGT");
        assert_eq!(record.quality_scores(), b"@DLS");
        Ok(())
    }

    #[test]
    fn test_read_record_quality_contains_plus() -> io::Result<()> {
        let data = b"@r0\nACGT\n+\n+DLS\n";
        let mut reader = &data[..];
        let mut record = Record::default();

        read_record(&mut reader, &mut record)?;
        assert_eq!(record.sequence(), b"ACGT");
        assert_eq!(record.quality_scores(), b"+DLS");
        Ok(())
    }

    #[test]
    fn test_read_record_no_trailing_newline() -> io::Result<()> {
        let data = b"@r0\nACGT\n+\nNDLS";
        let mut reader = &data[..];
        let mut record = Record::default();

        read_record(&mut reader, &mut record)?;
        assert_eq!(record.name(), &b"r0"[..]);
        assert_eq!(record.sequence(), b"ACGT");
        assert_eq!(record.quality_scores(), b"NDLS");
        Ok(())
    }

    #[test]
    fn test_read_record_fallback_small_buffer() -> io::Result<()> {
        use std::io::BufReader;

        let data = b"@r0\nACGTACGT\n+\nNDLSNDLS\n";
        let mut reader = BufReader::with_capacity(8, &data[..]);
        let mut record = Record::default();

        read_record(&mut reader, &mut record)?;
        assert_eq!(record.name(), &b"r0"[..]);
        assert_eq!(record.sequence(), b"ACGTACGT");
        assert_eq!(record.quality_scores(), b"NDLSNDLS");
        Ok(())
    }
}
