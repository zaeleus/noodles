use std::io::{self, BufRead};

use crate::Line;

pub(super) fn read_line<R>(reader: &mut R, line: &mut Line) -> io::Result<usize>
where
    R: BufRead,
{
    let dst = &mut line.0;
    dst.clear();
    super::read_line(reader, dst)
}
