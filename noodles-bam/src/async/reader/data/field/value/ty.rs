use std::convert::TryFrom;

use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::record::data::field::value::Type;

pub async fn read_type<R>(reader: &mut R) -> io::Result<Type>
where
    R: AsyncRead + Unpin,
{
    reader
        .read_u8()
        .await
        .and_then(|n| Type::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_type() -> io::Result<()> {
        let data = [b'i'];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader).await?, Type::Int32);

        let data = [b'n'];
        let mut reader = &data[..];
        assert!(matches!(
            read_type(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
