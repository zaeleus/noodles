//! Async BGZF I/O.

mod block_codec;
pub mod reader;
pub mod writer;

pub use self::{reader::Reader, writer::Writer};

use self::block_codec::BlockCodec;

#[cfg(test)]
mod tests {
    use tokio::io::{self, AsyncReadExt, AsyncWriteExt};

    use super::*;

    #[tokio::test]
    async fn test_self() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());
        writer.write_all(b"noodles").await?;
        writer.shutdown().await?;

        let data = writer.into_inner();
        let mut reader = Reader::new(&data[..]);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).await?;

        assert_eq!(buf, b"noodles");

        Ok(())
    }
}
