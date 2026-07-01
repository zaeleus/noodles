mod reference_sequence_names;

use noodles_csi::binning_index::index::{Header, header::Format};
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

    write_end_position_index(
        writer,
        header.format(),
        header.start_position_index(),
        header.end_position_index(),
    )
    .await?;

    let meta = i32::from(header.line_comment_prefix());
    writer.write_i32_le(meta).await?;

    let skip = i32::try_from(header.line_skip_count())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32_le(skip).await?;

    write_reference_sequence_names(writer, header.reference_sequence_names()).await?;

    Ok(())
}

async fn write_end_position_index<W>(
    writer: &mut W,
    format: Format,
    start_position_index: usize,
    end_position_index: Option<usize>,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    const SPECIALIZED_END_VALUE: i32 = 0;

    if matches!(format, Format::Sam | Format::Vcf) {
        if end_position_index.is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid end position index for format",
            ));
        } else {
            writer.write_i32_le(SPECIALIZED_END_VALUE).await?;
        }
    } else {
        let i = end_position_index.unwrap_or(start_position_index);
        let j = i.checked_add(1).expect("attempt to add with overflow");
        let n = i32::try_from(j).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        writer.write_i32_le(n).await?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_csi::binning_index::index::header::format::CoordinateSystem;

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
        async fn t(
            buf: &mut Vec<u8>,
            format: Format,
            start_position_index: usize,
            end_position_index: Option<usize>,
            expected: i32,
        ) -> io::Result<()> {
            buf.clear();
            write_end_position_index(buf, format, start_position_index, end_position_index).await?;
            assert_eq!(buf, &expected.to_le_bytes());
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, Format::Sam, 5, None, 0).await?;
        t(&mut buf, Format::Vcf, 5, None, 0).await?;
        t(&mut buf, Format::Generic(CoordinateSystem::Gff), 5, None, 6).await?;
        t(
            &mut buf,
            Format::Generic(CoordinateSystem::Gff),
            5,
            Some(8),
            9,
        )
        .await?;

        buf.clear();
        assert!(matches!(
            write_end_position_index(&mut buf, Format::Sam, 5, Some(8)).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        buf.clear();
        assert!(matches!(
            write_end_position_index(&mut buf, Format::Vcf, 5, Some(8)).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
