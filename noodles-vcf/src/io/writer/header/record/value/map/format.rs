mod number;

use std::io::{self, Write};

use self::number::write_number;
use super::{
    write_delimiter, write_description_field, write_key, write_other_fields, write_separator,
    write_type_field,
};
use crate::header::{
    record::value::{map::Format, Map},
    Number,
};

pub(crate) fn write_format<W>(writer: &mut W, format: &Map<Format>) -> io::Result<()>
where
    W: Write,
{
    write_number_field(writer, format.number())?;
    write_type_field(writer, format.ty())?;
    write_description_field(writer, format.description())?;
    write_other_fields(writer, format.other_fields())?;
    Ok(())
}

fn write_number_field<W>(writer: &mut W, number: Number) -> io::Result<()>
where
    W: Write,
{
    use crate::header::record::value::map::tag::NUMBER;

    write_delimiter(writer)?;
    write_key(writer, NUMBER)?;
    write_separator(writer)?;
    write_number(writer, number)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record::samples::keys::key;

    #[test]
    fn test_write_format() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        buf.clear();
        let map = Map::<Format>::from(key::GENOTYPE);
        write_format(&mut buf, &map)?;
        assert_eq!(buf, br#",Number=1,Type=String,Description="Genotype""#);

        buf.clear();
        let mut map = Map::<Format>::from(key::GENOTYPE);
        map.other_fields_mut()
            .insert("noodles".parse()?, String::from("vcf"));
        write_format(&mut buf, &map)?;
        assert_eq!(
            buf,
            br#",Number=1,Type=String,Description="Genotype",noodles="vcf""#
        );

        Ok(())
    }
}
