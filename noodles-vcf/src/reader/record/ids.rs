mod id;

use std::io;

use self::id::parse_id;
use crate::record::Ids;

pub(super) fn parse_ids(s: &str, ids: &mut Ids) -> io::Result<()> {
    const DELIMITER: char = ';';

    if s.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "empty IDs"));
    }

    for raw_id in s.split(DELIMITER) {
        let id = parse_id(raw_id).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        if !ids.insert(id) {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "duplicate ID"));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_ids() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::ids::Id;

        let id0: Id = "nd0".parse()?;
        let id1: Id = "nd1".parse()?;

        let mut ids = Ids::default();

        ids.clear();
        parse_ids("nd0", &mut ids)?;
        let expected = [id0.clone()].into_iter().collect();
        assert_eq!(ids, expected);

        ids.clear();
        parse_ids("nd0;nd1", &mut ids)?;
        let expected = [id0, id1].into_iter().collect();
        assert_eq!(ids, expected);

        ids.clear();
        assert!(matches!(
            parse_ids("", &mut ids),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        ids.clear();
        assert!(matches!(
            parse_ids("nd0;nd0", &mut ids),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        ids.clear();
        assert!(matches!(
            parse_ids("nd 0", &mut ids),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        ids.clear();
        assert!(matches!(
            parse_ids(";nd0", &mut ids),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        ids.clear();
        assert!(matches!(
            parse_ids("nd0;;nd1", &mut ids),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        ids.clear();
        assert!(matches!(
            parse_ids("nd0;", &mut ids),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
