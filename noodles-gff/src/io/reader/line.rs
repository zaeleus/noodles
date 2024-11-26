use std::io::{self, BufRead};

use crate::Line;

pub(super) fn read_line<R>(reader: &mut R, line: &mut Line) -> io::Result<usize>
where
    R: BufRead,
{
    let buf = &mut line.0;

    loop {
        buf.clear();

        let n = super::read_line(reader, buf)?;

        if n == 0 || !is_blank(buf) {
            return Ok(n);
        }
    }
}

pub(crate) fn is_blank(s: &str) -> bool {
    s.chars().all(char::is_whitespace)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_line() -> io::Result<()> {
        const DATA: &[u8] = b"\n#comment\n\t\n";

        let mut line = Line::default();
        let mut lines: Vec<String> = Vec::new();

        let mut src = DATA;

        while read_line(&mut src, &mut line)? != 0 {
            lines.push(line.as_ref().into());
        }

        assert_eq!(lines, [String::from("#comment")]);

        Ok(())
    }
}
