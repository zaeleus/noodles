mod record;

use std::io::{self, Write};

use self::record::{
    write_alternative_allele, write_contig, write_file_format, write_filter, write_format,
    write_info, write_other,
};
use crate::{header::SampleNames, Header};

pub(super) fn write_header<W>(writer: &mut W, header: &Header) -> io::Result<()>
where
    W: Write,
{
    write_file_format(writer, header.file_format())?;

    for (id, info) in header.infos() {
        write_info(writer, id, info)?;
    }

    for (id, filter) in header.filters() {
        write_filter(writer, id, filter)?;
    }

    for (id, format) in header.formats() {
        write_format(writer, id, format)?;
    }

    for (id, alternative_allele) in header.alternative_alleles() {
        write_alternative_allele(writer, id, alternative_allele)?;
    }

    for (id, contig) in header.contigs() {
        write_contig(writer, id, contig)?;
    }

    for (key, collection) in header.other_records() {
        write_other(writer, key, collection)?;
    }

    write_column_names(writer, header.sample_names())?;

    Ok(())
}

fn write_column_names<W>(writer: &mut W, sample_names: &SampleNames) -> io::Result<()>
where
    W: Write,
{
    fn write_delimiter<W>(writer: &mut W) -> io::Result<()>
    where
        W: Write,
    {
        const DELIMITER: u8 = b'\t';
        writer.write_all(&[DELIMITER])
    }

    const REQUIRED_FIELD_NAMES: [&str; 8] = [
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
    ];
    const FORMAT_FIELD_NAME: &str = "FORMAT";

    writer.write_all(REQUIRED_FIELD_NAMES[0].as_bytes())?;

    for name in REQUIRED_FIELD_NAMES.iter().skip(1) {
        write_delimiter(writer)?;
        writer.write_all(name.as_bytes())?;
    }

    if !sample_names.is_empty() {
        write_delimiter(writer)?;
        writer.write_all(FORMAT_FIELD_NAME.as_bytes())?;

        for sample_name in sample_names {
            write_delimiter(writer)?;
            writer.write_all(sample_name.as_bytes())?;
        }
    }

    write_newline(writer)?;

    Ok(())
}

fn write_newline<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const LINE_FEED: u8 = b'\n';
    writer.write_all(&[LINE_FEED])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_header() -> io::Result<()> {
        let mut buf = Vec::new();

        let header = Header::default();
        write_header(&mut buf, &header)?;

        let expected = b"##fileformat=VCFv4.4
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_column_names() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        let sample_names = SampleNames::new();
        write_column_names(&mut buf, &sample_names)?;
        assert_eq!(buf, b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

        buf.clear();
        let sample_names = [String::from("sample0")].into_iter().collect();
        write_column_names(&mut buf, &sample_names)?;
        assert_eq!(
            buf,
            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample0\n"
        );

        Ok(())
    }
}
