use tokio::io::{self, AsyncBufRead, AsyncBufReadExt};

use crate::Record;

pub(super) async fn read_record<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    record: &mut Record,
) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    const LINE_FEED: u8 = b'\n';

    buf.clear();

    if reader.read_until(LINE_FEED, buf).await? == 0 {
        return Ok(0);
    }

    let mut src = &buf[..];
    crate::io::reader::read_record(&mut src, record)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::fields::Bounds;

    #[tokio::test]
    async fn test_read_record() -> io::Result<()> {
        let mut buf = Vec::new();

        let mut src = &b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n"[..];
        let mut record = Record::default();
        read_record(&mut src, &mut buf, &mut record).await?;
        assert_eq!(record.fields().buf, b"*4*0255**00**");
        assert_eq!(record.fields().bounds, Bounds::default());

        let mut src = &b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\r\n"[..];
        let mut record = Record::default();
        read_record(&mut src, &mut buf, &mut record).await?;
        assert_eq!(record.fields().buf, b"*4*0255**00**");
        assert_eq!(record.fields().bounds, Bounds::default());

        let mut src = &b"\n"[..];
        assert!(matches!(
            read_record(&mut src, &mut buf, &mut record).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
