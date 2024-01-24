use std::{ffi::CString, num::NonZeroUsize};

use noodles_sam as sam;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

pub(super) async fn write_header<W>(writer: &mut W, header: &sam::Header) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    write_raw_header(writer, header).await?;
    write_reference_sequences(writer, header.reference_sequences()).await?;
    Ok(())
}

async fn write_raw_header<W>(writer: &mut W, header: &sam::Header) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::MAGIC_NUMBER;

    writer.write_all(MAGIC_NUMBER).await?;

    let text = serialize_header(header)?;
    let l_text =
        u32::try_from(text.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(l_text).await?;

    writer.write_all(&text).await?;

    Ok(())
}

fn serialize_header(header: &sam::Header) -> io::Result<Vec<u8>> {
    let mut writer = sam::io::Writer::new(Vec::new());
    writer.write_header(header)?;
    Ok(writer.into_inner())
}

async fn write_reference_sequences<W>(
    writer: &mut W,
    reference_sequences: &sam::header::ReferenceSequences,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let n_ref = u32::try_from(reference_sequences.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(n_ref).await?;

    for (name, reference_sequence) in reference_sequences {
        write_reference_sequence(writer, name, reference_sequence.length()).await?;
    }

    Ok(())
}

async fn write_reference_sequence<W>(
    writer: &mut W,
    reference_sequence_name: &[u8],
    length: NonZeroUsize,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let c_name = CString::new(reference_sequence_name)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let name = c_name.as_bytes_with_nul();

    let l_name =
        u32::try_from(name.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(l_name).await?;
    writer.write_all(name).await?;

    let l_ref = u32::try_from(usize::from(length))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(l_ref).await?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_reference_sequence() -> io::Result<()> {
        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        let mut buf = Vec::new();
        write_reference_sequence(&mut buf, b"sq0", SQ0_LN).await?;

        let expected = [
            0x04, 0x00, 0x00, 0x00, // l_name = 4
            0x73, 0x71, 0x30, 0x00, // name = b"sq0\x00"
            0x08, 0x00, 0x00, 0x00, // l_ref = 8
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
