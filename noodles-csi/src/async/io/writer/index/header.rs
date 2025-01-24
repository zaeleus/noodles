use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::binning_index::index::Header;

pub(super) async fn write_header<W>(
    writer: &mut W,
    min_shift: u8,
    depth: u8,
    header: Option<&Header>,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    writer.write_i32_le(i32::from(min_shift)).await?;
    writer.write_i32_le(i32::from(depth)).await?;
    write_aux(writer, header).await?;
    Ok(())
}

pub(super) async fn write_aux<W>(writer: &mut W, header: Option<&Header>) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::io::writer::index::header::write_header as write_tabix_header;

    let mut aux = Vec::new();

    if let Some(hdr) = header {
        write_tabix_header(&mut aux, hdr)?;
    }

    let l_aux =
        i32::try_from(aux.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32_le(l_aux).await?;

    writer.write_all(&aux).await?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_header() -> io::Result<()> {
        let mut buf = Vec::new();
        write_header(&mut buf, 14, 5, None).await?;

        let expected = [
            0x0e, 0x00, 0x00, 0x00, // min_shift = 14
            0x05, 0x00, 0x00, 0x00, // depth = 5
            0x00, 0x00, 0x00, 0x00, // l_aux = 0
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
