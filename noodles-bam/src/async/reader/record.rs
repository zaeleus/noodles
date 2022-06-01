//! Async BAM record field readers.

use bytes::BytesMut;
use noodles_sam::alignment::Record;
use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(super) async fn read_record<R>(
    reader: &mut R,
    buf: &mut BytesMut,
    record: &mut Record,
) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    use crate::reader::record::decode_record;

    let block_size = match read_block_size(reader).await? {
        0 => return Ok(0),
        n => n,
    };

    buf.resize(block_size, Default::default());
    reader.read_exact(buf).await?;

    decode_record(buf, record)?;

    Ok(block_size)
}

pub(super) async fn read_block_size<R>(reader: &mut R) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    match reader.read_u32_le().await {
        Ok(bs) => usize::try_from(bs).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(0),
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_block_size() -> io::Result<()> {
        let data = [0x08, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_block_size(&mut reader).await?, 8);

        let data = [];
        let mut reader = &data[..];
        assert_eq!(read_block_size(&mut reader).await?, 0);

        Ok(())
    }

    #[tokio::test]
    async fn test_read_record() -> io::Result<()> {
        let data = [
            0x22, 0x00, 0x00, 0x00, // block_size = 34
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x02, // l_read_name = 2
            0xff, // mapq = 255
            0x48, 0x12, // bin = 4680
            0x00, 0x00, // n_cigar_op = 0
            0x04, 0x00, // flag = 4
            0x00, 0x00, 0x00, 0x00, // l_seq = 0
            0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
            0xff, 0xff, 0xff, 0xff, // next_pos = -1
            0x00, 0x00, 0x00, 0x00, // tlen = 0
            0x2a, 0x00, // read_name = "*\x00"
        ];

        let mut reader = &data[..];
        let mut buf = BytesMut::new();
        let mut record = Record::default();
        let block_size = read_record(&mut reader, &mut buf, &mut record).await?;

        assert_eq!(block_size, 34);
        assert_eq!(record, Record::default());

        Ok(())
    }
}
