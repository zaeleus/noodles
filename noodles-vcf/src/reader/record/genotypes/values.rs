use std::io;

use crate::{
    header::record::value::{map::Format, Map},
    reader::record::MISSING,
    record::genotypes::{sample::Value, Keys},
    Header,
};

pub(super) fn parse_values(
    header: &Header,
    keys: &Keys,
    s: &str,
    values: &mut Vec<Option<Value>>,
) -> io::Result<()> {
    const DELIMITER: char = ':';

    if s.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "empty values"));
    } else if s == MISSING {
        return Ok(());
    }

    let mut raw_values = s.split(DELIMITER);

    for (key, raw_value) in keys.iter().zip(&mut raw_values) {
        let value = if let Some(format) = header.formats().get(key) {
            parse_value(format, raw_value)?
        } else {
            let format = Map::<Format>::from(key);
            parse_value(&format, raw_value)?
        };

        values.push(value);
    }

    if raw_values.next().is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unexpected genotype value",
        ));
    }

    Ok(())
}

fn parse_value(format: &Map<Format>, s: &str) -> io::Result<Option<Value>> {
    match s {
        MISSING => Ok(None),
        _ => Value::from_str_format(s, format)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_values() -> Result<(), Box<dyn std::error::Error>> {
        use crate::header::format::key;

        let header = Header::default();
        let mut values = Vec::new();

        let keys = Keys::try_from(vec![key::GENOTYPE])?;
        values.clear();
        parse_values(&header, &keys, ".", &mut values)?;
        assert!(values.is_empty());

        let keys = Keys::try_from(vec![key::GENOTYPE])?;
        values.clear();
        parse_values(&header, &keys, "0|0", &mut values)?;
        assert_eq!(values, vec![Some(Value::String(String::from("0|0")))]);

        let keys = Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?;
        values.clear();
        parse_values(&header, &keys, "0|0:13", &mut values)?;
        assert_eq!(
            values,
            vec![
                Some(Value::String(String::from("0|0"))),
                Some(Value::Integer(13)),
            ]
        );

        let keys = Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?;
        values.clear();
        parse_values(&header, &keys, "0|0:.", &mut values)?;
        assert_eq!(values, vec![Some(Value::String(String::from("0|0"))), None]);

        let keys = Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?;
        values.clear();
        parse_values(&header, &keys, "0|0", &mut values)?;
        assert_eq!(values, vec![Some(Value::String(String::from("0|0")))]);

        let keys = Keys::try_from(vec![key::GENOTYPE])?;
        values.clear();
        assert!(matches!(
            parse_values(&header, &keys, "", &mut values),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        let keys = Keys::try_from(vec![key::GENOTYPE])?;
        values.clear();
        assert!(matches!(
            parse_values(&header, &keys, "0|0:13", &mut values),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
