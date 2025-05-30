use std::io::{self, Write};

use super::{write_description_field, write_other_fields};
use crate::header::record::value::{Map, map::Filter};

pub(crate) fn write_filter<W>(writer: &mut W, filter: &Map<Filter>) -> io::Result<()>
where
    W: Write,
{
    write_description_field(writer, filter.description())?;
    write_other_fields(writer, filter.other_fields())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_filter() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        buf.clear();
        let map = Map::<Filter>::new("All filters passed");
        write_filter(&mut buf, &map)?;
        assert_eq!(buf, br#",Description="All filters passed""#);

        buf.clear();
        let map = Map::<Filter>::builder()
            .set_description("All filters passed")
            .insert("noodles".parse()?, "vcf")
            .build()?;
        write_filter(&mut buf, &map)?;
        assert_eq!(buf, br#",Description="All filters passed",noodles="vcf""#);

        Ok(())
    }
}
