use noodles_csi::binning_index::index::{header::ReferenceSequenceNames, Header};
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

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

    let col_end = header.end_position_index().map_or(Ok(0), |mut i| {
        i = i.checked_add(1).expect("attempt to add with overflow");
        i32::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    })?;
    writer.write_i32_le(col_end).await?;

    let meta = i32::from(header.line_comment_prefix());
    writer.write_i32_le(meta).await?;

    let skip = i32::try_from(header.line_skip_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32_le(skip).await?;

    write_reference_sequence_names(writer, header.reference_sequence_names()).await?;

    Ok(())
}

async fn write_reference_sequence_names<W>(
    writer: &mut W,
    reference_sequence_names: &ReferenceSequenceNames,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    const NUL: u8 = b'\x00';

    // Add 1 for each trailing NUL.
    let len: usize = reference_sequence_names.iter().map(|n| n.len() + 1).sum();
    let l_nm = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32_le(l_nm).await?;

    for reference_sequence_name in reference_sequence_names {
        writer.write_all(reference_sequence_name).await?;
        writer.write_u8(NUL).await?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use bstr::BString;

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
    async fn test_write_reference_sequence_names() -> io::Result<()> {
        let reference_sequence_names = [BString::from("sq0"), BString::from("sq1")]
            .into_iter()
            .collect();

        let mut buf = Vec::new();
        write_reference_sequence_names(&mut buf, &reference_sequence_names).await?;

        let expected = [
            0x08, 0x00, 0x00, 0x00, // l_nm = 8
            0x73, 0x71, 0x30, 0x00, // names[0] = b"sq0\x00"
            0x73, 0x71, 0x31, 0x00, // names[1] = b"sq1\x00"
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
