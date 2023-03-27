use std::io;

use crate::{
    record::{info::field, Info},
    Header,
};

pub(super) fn parse_info(header: &Header, s: &str, info: &mut Info) -> io::Result<()> {
    const DELIMITER: char = ';';

    if s.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "empty info"));
    }

    for raw_field in s.split(DELIMITER) {
        let (key, value) = field::parse(raw_field, header.infos())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        if info.insert(key, value).is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "duplicate info key",
            ));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_info() -> io::Result<()> {
        use crate::{header::info::key, record::info::field::Value};

        let header = Header::default();

        let mut info = Info::default();

        info.clear();
        parse_info(&header, "NS=2", &mut info)?;
        let expected = [(key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2)))]
            .into_iter()
            .collect();
        assert_eq!(info, expected);

        info.clear();
        parse_info(&header, "NS=2;AA=T", &mut info)?;
        let expected = [
            (key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2))),
            (
                key::ANCESTRAL_ALLELE,
                Some(Value::String(String::from("T"))),
            ),
        ]
        .into_iter()
        .collect();
        assert_eq!(info, expected);

        info.clear();
        assert!(matches!(
            parse_info(&header, ".", &mut info),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        info.clear();
        assert!(matches!(
            parse_info(&header, "NS=ndls", &mut info),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
