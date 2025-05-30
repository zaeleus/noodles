use std::io::{self, Write};

use super::{write_delimiter, write_other_fields, write_value_field};
use crate::header::record::value::{
    Map,
    map::{Contig, contig::tag},
};

pub(crate) fn write_contig<W>(writer: &mut W, contig: &Map<Contig>) -> io::Result<()>
where
    W: Write,
{
    if let Some(length) = contig.length() {
        write_delimiter(writer)?;
        write_value_field(writer, tag::LENGTH, length.to_string())?;
    }

    if let Some(md5) = contig.md5() {
        write_delimiter(writer)?;
        write_value_field(writer, tag::MD5, md5)?;
    }

    if let Some(url) = contig.url() {
        write_delimiter(writer)?;
        write_value_field(writer, tag::URL, url)?;
    }

    write_other_fields(writer, contig.other_fields())?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_contig() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        buf.clear();
        let map = Map::<Contig>::new();
        write_contig(&mut buf, &map)?;
        assert!(buf.is_empty());

        buf.clear();
        let map = Map::<Contig>::builder()
            .set_length(8)
            .set_md5("d7eba311421bbc9d3ada44709dd61534")
            .set_url("https://example.com/reference.fa")
            .insert("noodles".parse()?, "vcf")
            .build()?;
        write_contig(&mut buf, &map)?;
        assert_eq!(
            buf,
            br#",length=8,md5=d7eba311421bbc9d3ada44709dd61534,URL=https://example.com/reference.fa,noodles="vcf""#
        );

        Ok(())
    }
}
