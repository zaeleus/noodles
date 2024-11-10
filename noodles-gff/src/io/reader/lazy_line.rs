use std::io::{self, BufRead};

use super::read_line;
use crate::lazy;

pub(super) fn read_lazy_line<R>(reader: &mut R, line: &mut lazy::Line) -> io::Result<usize>
where
    R: BufRead,
{
    let buf = &mut line.0;
    buf.clear();
    read_line(reader, buf)
}
