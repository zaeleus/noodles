use std::io::{self, Write};

use super::{write_description_field, write_other_fields};
use crate::header::record::value::{Map, map::AlternativeAllele};

pub(crate) fn write_alternative_allele<W>(
    writer: &mut W,
    alternative_allele: &Map<AlternativeAllele>,
) -> io::Result<()>
where
    W: Write,
{
    write_description_field(writer, alternative_allele.description())?;
    write_other_fields(writer, alternative_allele.other_fields())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_alternative_allele() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        buf.clear();
        let map = Map::<AlternativeAllele>::new("Deletion");
        write_alternative_allele(&mut buf, &map)?;
        assert_eq!(buf, br#",Description="Deletion""#);

        buf.clear();
        let map = Map::<AlternativeAllele>::builder()
            .set_description("Deletion")
            .insert("noodles".parse()?, "vcf")
            .build()?;
        write_alternative_allele(&mut buf, &map)?;
        assert_eq!(buf, br#",Description="Deletion",noodles="vcf""#);

        Ok(())
    }
}
