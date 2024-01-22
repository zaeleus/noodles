use std::{
    io::{self, BufRead},
    str,
};

use super::read_line;
use crate::lazy;

pub(super) fn read_lazy_record<R>(reader: &mut R, record: &mut lazy::Record) -> io::Result<usize>
where
    R: BufRead,
{
    record.buf.clear();

    let mut len = 0;

    len += read_field(reader, &mut record.buf)?;
    record.bounds.chromosome_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.position_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.ids_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.reference_bases_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.alternate_bases_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.quality_score_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.filters_end = record.buf.len();

    len += read_field(reader, &mut record.buf)?;
    record.bounds.info_end = record.buf.len();

    len += read_line(reader, &mut record.buf)?;

    Ok(len)
}

fn read_field<R>(reader: &mut R, dst: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    const DELIMITER: u8 = b'\t';

    let mut is_delimiter = false;
    let mut len = 0;

    loop {
        let src = reader.fill_buf()?;

        if is_delimiter || src.is_empty() {
            break;
        }

        let (buf, n) = match src.iter().position(|&b| b == DELIMITER) {
            Some(i) => {
                is_delimiter = true;
                (&src[..i], i + 1)
            }
            None => (src, src.len()),
        };

        let s = str::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        dst.push_str(s);

        len += n;

        reader.consume(n);
    }

    Ok(len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_lazy_record() -> io::Result<()> {
        let mut src = &b"sq0\t1\t.\tA\t.\t.\tPASS\t."[..];

        let mut record = lazy::Record::default();
        read_lazy_record(&mut src, &mut record)?;

        assert_eq!(record.buf, "sq01.A..PASS.");

        assert_eq!(record.bounds.chromosome_end, 3);
        assert_eq!(record.bounds.position_end, 4);
        assert_eq!(record.bounds.ids_end, 5);
        assert_eq!(record.bounds.reference_bases_end, 6);
        assert_eq!(record.bounds.alternate_bases_end, 7);
        assert_eq!(record.bounds.quality_score_end, 8);
        assert_eq!(record.bounds.filters_end, 12);
        assert_eq!(record.bounds.info_end, 13);

        Ok(())
    }
}
