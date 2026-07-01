mod reference_sequence_names;

use noodles_csi::binning_index::index::Header;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::reference_sequence_names::write_reference_sequence_names;

pub(super) async fn write_header<W>(writer: &mut W, header: &Header) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let format = i32::from(header.format());
    writer.write_i32_le(format).await?;

    let reference_sequence_name_index = header
        .reference_sequence_name_index()
        .checked_add(1)
        .expect("attempt to add with overflow");
    let col_seq = i32::try_from(reference_sequence_name_index)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32_le(col_seq).await?;

    let start_position_index = header
        .start_position_index()
        .checked_add(1)
        .expect("attempt to add with overflow");
    let col_beg = i32::try_from(start_position_index)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32_le(col_beg).await?;

    write_end_position_index(writer, header.end_position_index()).await?;

    let meta = i32::from(header.line_comment_prefix());
    writer.write_i32_le(meta).await?;

    let skip = i32::try_from(header.line_skip_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32_le(skip).await?;

    write_reference_sequence_names(writer, header.reference_sequence_names()).await?;

    Ok(())
}

async fn write_end_position_index<W>(writer: &mut W, i: Option<usize>) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let n = i.map_or(Ok(0), |mut j| {
        j = j.checked_add(1).expect("attempt to add with overflow");
        i32::try_from(j).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    })?;

    writer.write_i32_le(n).await
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_header() -> io::Result<()> {
        use noodles_csi::binning_index::index::header;

        let header = header::Builder::gff().build();

        let mut buf = Vec::new();
        write_header(&mut buf, &header).await?;

        let expected = [
            0x00, 0x00, 0x00, 0x00, // format = Generic(GFF)
            0x01, 0x00, 0x00, 0x00, // col_seq = 1
            0x04, 0x00, 0x00, 0x00, // col_beg = 4
            0x05, 0x00, 0x00, 0x00, // col_end = 5
            0x23, 0x00, 0x00, 0x00, // meta = '#'
            0x00, 0x00, 0x00, 0x00, // skip = 0
            0x00, 0x00, 0x00, 0x00, // l_nm = 0
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[tokio::test]
    async fn test_write_end_position_index() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_end_position_index(&mut buf, None).await?;
        assert_eq!(buf, 0i32.to_le_bytes());

        buf.clear();
        write_end_position_index(&mut buf, Some(0)).await?;
        assert_eq!(buf, 1i32.to_le_bytes());

        Ok(())
    }
}
