use std::io;

use crate::record::ReadName;

const MAX_LENGTH: usize = 254;

pub(crate) fn parse_read_name(src: &[u8]) -> io::Result<ReadName> {
    if src.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "empty read name",
        ));
    } else if src.len() > MAX_LENGTH {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "read name too long",
        ));
    } else if !is_valid_name(src) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid read name",
        ));
    }

    ReadName::try_new(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn is_valid_name_char(b: u8) -> bool {
    matches!(b, b'!'..=b'?' | b'A'..=b'~')
}

fn is_valid_name(s: &[u8]) -> bool {
    s.iter().copied().all(is_valid_name_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_name() -> Result<(), Box<dyn std::error::Error>> {
        assert_eq!(parse_read_name(b"r0")?, ReadName::try_new("r0")?);

        assert!(matches!(
            parse_read_name(b""),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        assert!(matches!(
            parse_read_name(b"r 0"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        assert!(matches!(
            parse_read_name(b"@r0"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        let src = vec![b'n'; MAX_LENGTH + 1];
        assert!(matches!(
            parse_read_name(&src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
