use std::io::{self, Write};

use super::{write_description_field, write_number_field, write_other_fields, write_type_field};
use crate::header::record::value::{map::Format, Map};

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::genotypes::keys::key;

    #[test]
    fn test_write_format() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        buf.clear();
        let map = Map::<Format>::from(&key::GENOTYPE);
        write_format(&mut buf, &map)?;
        assert_eq!(buf, br#",Number=1,Type=String,Description="Genotype""#);

        buf.clear();
        let mut map = Map::<Format>::from(&key::GENOTYPE);
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
