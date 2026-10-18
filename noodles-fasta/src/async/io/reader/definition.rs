use tokio::io::{self, AsyncBufRead};

use super::read_line;
use crate::{
    io::reader::definition::{is_comment, parse_definition},
    record::Definition,
};

pub(super) async fn read_definition<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    definition: &mut Definition,
) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    let mut len = 0;

    loop {
        buf.clear();

        match read_line(reader, buf).await? {
            0 => return Ok(0),
            n if is_comment(buf) => len += n,
            n => {
                len += n;
                break;
            }
        }
    }

    let (name, description) = parse_definition(buf)?;

    let description = if description.is_empty() {
        None
    } else {
        Some(description.into())
    };

    *definition = Definition::new(name, description);

    Ok(len)
}
