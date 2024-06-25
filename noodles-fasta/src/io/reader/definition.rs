use std::io::{self, BufRead};

use super::read_line;

pub(super) fn read_definition<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    read_line(reader, buf)
}
