use tokio::io::{self, AsyncBufRead, AsyncBufReadExt};

use crate::Record;

pub(super) async fn read_record<R>(
    reader: &mut R,
    buf: &mut String,
    record: &mut Record,
) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    buf.clear();

    if reader.read_line(buf).await? == 0 {
        return Ok(0);
    }

    let mut buf = buf.as_bytes();
    crate::io::reader::record::read_record(&mut buf, record)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_record() -> io::Result<()> {
        let mut src = &b"sq0\t1\t.\tA\t.\t.\tPASS\t.\n"[..];
        let mut buf = String::new();

        let mut record = Record::default();
        read_record(&mut src, &mut buf, &mut record).await?;

        assert_eq!(record.fields().buf, "sq01.A..PASS.");

        assert_eq!(record.fields().bounds.reference_sequence_name_end, 3);
        assert_eq!(record.fields().bounds.variant_start_end, 4);
        assert_eq!(record.fields().bounds.ids_end, 5);
        assert_eq!(record.fields().bounds.reference_bases_end, 6);
        assert_eq!(record.fields().bounds.alternate_bases_end, 7);
        assert_eq!(record.fields().bounds.quality_score_end, 8);
        assert_eq!(record.fields().bounds.filters_end, 12);
        assert_eq!(record.fields().bounds.info_end, 13);

        let mut src = &b"\n"[..];
        assert!(matches!(
            read_record(&mut src, &mut buf, &mut record).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
