use std::io;

use crate::record::Position;

pub(super) fn parse_position(s: &str) -> io::Result<Position> {
    let n: usize = s
        .parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(Position::from(n))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_position() -> io::Result<()> {
        assert_eq!(parse_position("0")?, Position::from(0));
        assert_eq!(parse_position("8")?, Position::from(8));

        assert!(matches!(
            parse_position(""),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));
        assert!(matches!(
            parse_position("ndls"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));
        assert!(matches!(
            parse_position("-1"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
