use noodles_csi::binning_index::index::header::ReferenceSequenceNames;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

pub(super) async fn write_reference_sequence_names<W>(
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
    async fn test_write_reference_sequence_names() -> io::Result<()> {
        let mut buf = Vec::new();

        let reference_sequence_names = ReferenceSequenceNames::default();
        buf.clear();
        write_reference_sequence_names(&mut buf, &reference_sequence_names).await?;
        assert_eq!(buf, [0x00, 0x00, 0x00, 0x00]); // l_nm = 0

        let reference_sequence_names = [BString::from("sq0")].into_iter().collect();
        buf.clear();
        write_reference_sequence_names(&mut buf, &reference_sequence_names).await?;
        assert_eq!(
            buf,
            [
                0x04, 0x00, 0x00, 0x00, // l_nm = 4
                b's', b'q', b'0', 0x00, // names[0] = "sq0"
            ]
        );

        let reference_sequence_names = [BString::from("sq0"), BString::from("sq1")]
            .into_iter()
            .collect();
        buf.clear();
        write_reference_sequence_names(&mut buf, &reference_sequence_names).await?;
        assert_eq!(
            buf,
            [
                0x08, 0x00, 0x00, 0x00, // l_nm = 8
                b's', b'q', b'0', 0x00, // names[0] = "sq0"
                b's', b'q', b'1', 0x00, // names[1] = "sq1"
            ]
        );

        Ok(())
    }
}
