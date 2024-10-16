use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::record::Definition;

pub(super) async fn write_definition<W>(writer: &mut W, definition: &Definition) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    write_prefix(writer).await?;
    write_name(writer, definition.name()).await?;

    if let Some(description) = definition.description() {
        write_separator(writer).await?;
        write_description(writer, description).await?;
    }

    Ok(())
}

async fn write_prefix<W>(writer: &mut W) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    const PREFIX: u8 = b'>';
    writer.write_all(&[PREFIX]).await
}

async fn write_name<W>(writer: &mut W, name: &[u8]) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    writer.write_all(name).await
}

async fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    const SEPARATOR: u8 = b' ';
    writer.write_all(&[SEPARATOR]).await
}

async fn write_description<W>(writer: &mut W, description: &[u8]) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    writer.write_all(description).await
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_definition() -> io::Result<()> {
        async fn t(buf: &mut Vec<u8>, definition: &Definition, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_definition(buf, definition).await?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &Definition::new("sq0", None), b">sq0").await?;
        t(
            &mut buf,
            &Definition::new("sq0", Some(Vec::from("LN:8"))),
            b">sq0 LN:8",
        )
        .await?;

        Ok(())
    }
}
