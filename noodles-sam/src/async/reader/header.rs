use tokio::io::{self, AsyncBufRead, AsyncBufReadExt};

pub(super) async fn read_header<R>(reader: &mut R) -> io::Result<String>
where
    R: AsyncBufRead + Unpin,
{
    const PREFIX: u8 = b'@';
    const LINE_FEED: u8 = b'\n';

    let mut buf = Vec::new();
    let mut is_eol = false;

    for i in 0.. {
        let src = reader.fill_buf().await?;

        if (i == 0 || is_eol) && src.first().map(|&b| b != PREFIX).unwrap_or(true) {
            break;
        }

        let (read_eol, len) = if let Some(i) = src.iter().position(|&b| b == LINE_FEED) {
            buf.extend(&src[..=i]);
            (true, i + 1)
        } else {
            buf.extend(src);
            (false, src.len())
        };

        is_eol = read_eol;

        reader.consume(len);
    }

    String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
