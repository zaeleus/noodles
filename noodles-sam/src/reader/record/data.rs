pub(crate) mod field;

use std::io;

use self::field::parse_field;
use crate::record::Data;

pub(crate) fn parse_data(mut src: &[u8], data: &mut Data) -> io::Result<()> {
    while !src.is_empty() {
        let (tag, value) = parse_field(&mut src)?;
        data.insert(tag, value);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_data() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::data::field::{Tag, Value};

        let mut data = Data::default();

        parse_data(b"", &mut data)?;
        assert!(data.is_empty());

        let nh = (Tag::AlignmentHitCount, Value::from(1u8));
        let co = (Tag::Comment, Value::String(String::from("ndls")));

        data.clear();
        parse_data(b"NH:i:1", &mut data)?;
        let expected = [nh.clone()].into_iter().collect();
        assert_eq!(data, expected);

        data.clear();
        parse_data(b"NH:i:1\tCO:Z:ndls", &mut data)?;
        let expected = [nh, co].into_iter().collect();
        assert_eq!(data, expected);

        Ok(())
    }
}
