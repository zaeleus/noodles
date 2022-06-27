pub(crate) mod field;

use std::io;

use self::field::parse_field;
use crate::record::Data;

pub(crate) fn parse_data(mut src: &[u8]) -> io::Result<Data> {
    let mut data = Data::default();

    while let Some(field) = parse_field(&mut src)? {
        data.insert(field);
    }

    Ok(data)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_data() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::data::{
            field::{Tag, Value},
            Field,
        };

        assert!(parse_data(b"")?.is_empty());

        let nh = Field::new(Tag::AlignmentHitCount, Value::from(1u8));
        let co = Field::new(Tag::Comment, Value::String(String::from("ndls")));

        let actual = parse_data(b"NH:i:1")?;
        let expected = Data::try_from(vec![nh.clone()])?;
        assert_eq!(actual, expected);

        let actual = parse_data(b"NH:i:1\tCO:Z:ndls")?;
        let expected = Data::try_from(vec![nh, co])?;
        assert_eq!(actual, expected);

        Ok(())
    }
}
