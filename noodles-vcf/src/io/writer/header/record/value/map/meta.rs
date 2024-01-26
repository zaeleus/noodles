use std::io::{self, Write};

use super::{write_delimiter, write_string_field, write_value_field};
use crate::header::record::value::{map::Other, Map};

pub(crate) fn write_meta<W>(writer: &mut W, meta: &Map<Other>) -> io::Result<()>
where
    W: Write,
{
    const NUMBER: &str = "Number";
    const TYPE: &str = "Type";
    const VALUES: &str = "Values";

    for (key, value) in meta.other_fields() {
        write_delimiter(writer)?;

        match key.as_ref() {
            NUMBER | TYPE | VALUES => write_value_field(writer, key, value)?,
            _ => write_string_field(writer, key, value)?,
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_meta() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        buf.clear();
        let map = Map::<Other>::builder()
            .insert("Type".parse()?, "String")
            .insert("Number".parse()?, ".")
            .insert("Values".parse()?, "[WholeGenome, Exome]")
            .insert("noodles".parse()?, "vcf")
            .build()?;
        write_meta(&mut buf, &map)?;
        assert_eq!(
            buf,
            br#",Type=String,Number=.,Values=[WholeGenome, Exome],noodles="vcf""#
        );

        Ok(())
    }
}
