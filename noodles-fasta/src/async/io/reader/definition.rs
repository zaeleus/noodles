use tokio::io::{self, AsyncBufRead};

use super::read_line;
use crate::{io::reader::definition::parse_definition, record::Definition};

pub(super) async fn read_definition<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    definition: &mut Definition,
) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    buf.clear();

    match read_line(reader, buf).await? {
        0 => Ok(0),
        n => {
            let (name, description) = parse_definition(buf)?;

            let description = if description.is_empty() {
                None
            } else {
                Some(description.into())
            };

            *definition = Definition::new(name, description);

            Ok(n)
        }
    }
}
