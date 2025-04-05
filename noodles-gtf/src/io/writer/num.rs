use std::io::{self, Write};

use lexical_core::FormattedSize;

pub fn write_usize<W>(writer: &mut W, n: usize) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; usize::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}
