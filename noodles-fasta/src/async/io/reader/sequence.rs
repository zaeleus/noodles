use tokio::io::{self, AsyncBufRead, AsyncBufReadExt};

pub(super) async fn read_sequence<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    use memchr::memchr;

    use crate::io::reader::DEFINITION_PREFIX;

    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    let mut n = 0;

    loop {
        let src = reader.fill_buf().await?;

        if src.first().map(|&b| b == DEFINITION_PREFIX).unwrap_or(true) {
            break;
        }

        let len = match memchr(LINE_FEED, src) {
            Some(i) => {
                let line = &src[..i];

                if line.ends_with(&[CARRIAGE_RETURN]) {
                    let end = line.len() - 1;
                    buf.extend_from_slice(&line[..end]);
                } else {
                    buf.extend_from_slice(line);
                }

                i + 1
            }
            None => {
                buf.extend(src);
                src.len()
            }
        };

        reader.consume(len);

        n += len;
    }

    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_sequence() -> io::Result<()> {
        async fn t(buf: &mut Vec<u8>, mut reader: &[u8], expected: &[u8]) -> io::Result<()> {
            buf.clear();
            read_sequence(&mut reader, buf).await?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, b"ACGT\n", b"ACGT").await?;
        t(&mut buf, b"ACGT\n>sq1\n", b"ACGT").await?;
        t(&mut buf, b"NNNN\nNNNN\nNN\n", b"NNNNNNNNNN").await?;

        Ok(())
    }
}
