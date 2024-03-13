use std::io::{self, Write};

use super::{write_description_field, write_number_field, write_other_fields, write_type_field};
use crate::header::record::value::{map::Info, Map};

pub(crate) fn write_info<W>(writer: &mut W, info: &Map<Info>) -> io::Result<()>
where
    W: Write,
{
    write_number_field(writer, info.number())?;
    write_type_field(writer, info.ty())?;
    write_description_field(writer, info.description())?;
    write_other_fields(writer, info.other_fields())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record::info::field::key;

    #[test]
    fn test_write_info() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        buf.clear();
        let map = Map::<Info>::from(key::SAMPLES_WITH_DATA_COUNT);
        write_info(&mut buf, &map)?;
        assert_eq!(
            buf,
            br#",Number=1,Type=Integer,Description="Number of samples with data""#
        );

        buf.clear();
        let mut map = Map::<Info>::from(key::SAMPLES_WITH_DATA_COUNT);
        map.other_fields_mut()
            .insert("noodles".parse()?, String::from("vcf"));
        write_info(&mut buf, &map)?;
        assert_eq!(
            buf,
            br#",Number=1,Type=Integer,Description="Number of samples with data",noodles="vcf""#
        );

        Ok(())
    }
}
