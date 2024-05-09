mod number;
mod ty;

use std::io::{self, Write};

use self::{number::write_number, ty::write_type};
use super::{
    write_delimiter, write_description_field, write_key, write_other_fields, write_separator,
};
use crate::header::record::value::{
    map::{
        format::{Number, Type},
        Format,
    },
    Map,
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

fn write_type_field<W>(writer: &mut W, ty: Type) -> io::Result<()>
where
    W: Write,
{
    use crate::header::record::value::map::tag::TYPE;

    write_delimiter(writer)?;
    write_key(writer, TYPE)?;
    write_separator(writer)?;
    write_type(writer, ty)?;

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
