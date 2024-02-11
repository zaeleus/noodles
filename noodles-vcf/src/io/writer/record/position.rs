use std::io::{self, Write};

use noodles_core::Position;

pub(super) fn write_position<W>(writer: &mut W, position: Option<Position>) -> io::Result<()>
where
    W: Write,
{
    let n = position.map(usize::from).unwrap_or_default();
    write!(writer, "{}", n)
}
