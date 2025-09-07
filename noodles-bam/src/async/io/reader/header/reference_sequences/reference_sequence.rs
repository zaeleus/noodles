use std::num::NonZero;

use bstr::BString;
use noodles_sam::header::record::value::{Map, map::ReferenceSequence};
use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::io::reader::bytes_with_nul_to_bstring;

pub(super) async fn read_reference_sequence<R>(
    reader: &mut R,
) -> io::Result<(BString, Map<ReferenceSequence>)>
where
    R: AsyncRead + Unpin,
{
    let name = read_name(reader).await?;

    let len = read_length(reader).await?;
    let reference_sequence = Map::<ReferenceSequence>::new(len);

    Ok((name, reference_sequence))
}

async fn read_name<R>(reader: &mut R) -> io::Result<BString>
where
    R: AsyncRead + Unpin,
{
    let l_name = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut c_name = vec![0; l_name];
    reader.read_exact(&mut c_name).await?;

    bytes_with_nul_to_bstring(&c_name)
}

async fn read_length<R>(reader: &mut R) -> io::Result<NonZero<usize>>
where
    R: AsyncRead + Unpin,
{
    reader.read_u32_le().await.and_then(|len| {
        usize::try_from(len)
            .and_then(NonZero::try_from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_reference_sequence() -> io::Result<()> {
        let data = [
            0x04, 0x00, 0x00, 0x00, // l_name = 4
            0x73, 0x71, 0x30, 0x00, // name = "sq0\x00"
            0x08, 0x00, 0x00, 0x00, // l_ref = 8
        ];

        let mut reader = &data[..];
        let actual = read_reference_sequence(&mut reader).await?;
        let expected = (
            BString::from("sq0"),
            Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
        );
        assert_eq!(actual, expected);

        Ok(())
    }
}
