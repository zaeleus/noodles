mod header;
mod magic_number;
mod reference_sequences;

use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::{
    header::write_aux, magic_number::write_magic_number,
    reference_sequences::write_reference_sequences,
};
use crate::{BinningIndex, Index, io::MAX_DEPTH};

pub(super) async fn write_index<W>(writer: &mut W, index: &Index) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    write_magic_number(writer).await?;

    write_min_shift(writer, index.min_shift()).await?;
    write_depth(writer, index.depth()).await?;

    write_aux(writer, index.header()).await?;
    write_reference_sequences(writer, index.depth(), index.reference_sequences()).await?;

    if let Some(n) = index.unplaced_unmapped_record_count() {
        writer.write_u64_le(n).await?;
    }

    Ok(())
}

async fn write_min_shift<W>(writer: &mut W, min_shift: u8) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let n = i32::from(min_shift);
    writer.write_i32_le(n).await
}

async fn write_depth<W>(writer: &mut W, depth: u8) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    if depth > MAX_DEPTH {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid depth"));
    }

    let n = i32::from(depth);
    writer.write_i32_le(n).await
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_min_shift() -> io::Result<()> {
        let mut buf = Vec::new();
        write_min_shift(&mut buf, 14).await?;
        assert_eq!(buf, [0x0e, 0x00, 0x00, 0x00]);
        Ok(())
    }

    #[tokio::test]
    async fn test_write_depth() -> io::Result<()> {
        async fn t(buf: &mut Vec<u8>, depth: u8, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_depth(buf, depth).await?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 0, &[0x00, 0x00, 0x00, 0x00]).await?;
        t(&mut buf, 5, &[0x05, 0x00, 0x00, 0x00]).await?;
        t(&mut buf, 9, &[0x09, 0x00, 0x00, 0x00]).await?;

        buf.clear();
        assert!(matches!(
            write_depth(&mut buf, 10).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
