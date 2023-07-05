use crate::record::data::field::value::base_modifications::group::Status;

pub(super) fn parse_status(src: &mut &[u8]) -> Option<Status> {
    if let Some((b, rest)) = src.split_first() {
        let status = match *b {
            b'.' => Status::Implicit,
            b'?' => Status::Explicit,
            _ => return None,
        };

        *src = rest;

        Some(status)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_status() {
        let mut src = &[][..];
        assert!(parse_status(&mut src).is_none());

        let mut src = &b"."[..];
        assert_eq!(parse_status(&mut src), Some(Status::Implicit));

        let mut src = &b"?"[..];
        assert_eq!(parse_status(&mut src), Some(Status::Explicit));

        let mut src = &b";"[..];
        assert!(parse_status(&mut src).is_none());
    }
}
