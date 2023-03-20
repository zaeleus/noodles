use std::io;

use crate::{reader::record::MISSING, record::genotypes::Keys, Header};

pub(super) fn parse_keys(header: &Header, s: &str, keys: &mut Keys) -> io::Result<()> {
    use crate::header::format::key;

    const DELIMITER: char = ':';

    if s.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "missing keys"));
    } else if s == MISSING {
        return Ok(());
    }

    let mut gt_position = None;

    for (i, raw_key) in s.split(DELIMITER).enumerate() {
        let key = match header.formats().get_full(raw_key) {
            Some((_, k, _)) => k.clone(),
            None => raw_key
                .parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        };

        if key == key::GENOTYPE {
            gt_position = Some(i);
        }

        if !keys.insert(key) {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "duplicate key"));
        }
    }

    if let Some(i) = gt_position {
        if i != 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid GT key position",
            ));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_keys() -> Result<(), Box<dyn std::error::Error>> {
        use crate::header::format::key;

        let header = Header::default();
        let mut keys = Keys::default();

        keys.clear();
        parse_keys(&header, ".", &mut keys)?;
        assert_eq!(keys, Keys::default());

        keys.clear();
        parse_keys(&header, "GT", &mut keys)?;
        let expected = Keys::try_from(vec![key::GENOTYPE])?;
        assert_eq!(keys, expected);

        keys.clear();
        parse_keys(&header, "GQ", &mut keys)?;
        let expected = Keys::try_from(vec![key::CONDITIONAL_GENOTYPE_QUALITY])?;
        assert_eq!(keys, expected);

        keys.clear();
        parse_keys(&header, "GT:GQ", &mut keys)?;
        let expected = Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?;
        assert_eq!(keys, expected);

        keys.clear();
        assert!(matches!(
            parse_keys(&header, "", &mut keys),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        keys.clear();
        assert!(matches!(
            parse_keys(&header, "GQ:GT", &mut keys),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        keys.clear();
        assert!(matches!(
            parse_keys(&header, "GT:GT", &mut keys),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
