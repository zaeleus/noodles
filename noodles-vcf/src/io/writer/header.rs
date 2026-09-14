mod record;

use std::{
    collections::HashSet,
    io::{self, Write},
};

use self::record::{
    write_alternative_allele, write_contig, write_file_format, write_filter, write_format,
    write_info, write_other,
};
use super::RecordOrder;
use crate::{
    Header,
    header::{RecordKey, SampleNames, record::key::Other},
};

pub(super) fn write_header<W>(
    writer: &mut W,
    header: &Header,
    record_order: RecordOrder,
) -> io::Result<()>
where
    W: Write,
{
    write_file_format(writer, header.file_format())?;

    match record_order {
        RecordOrder::Grouped => write_grouped_records(writer, header)?,
        RecordOrder::AsAdded => write_added_records(writer, header)?,
    }

    write_column_names(writer, header.sample_names())?;

    Ok(())
}

fn write_grouped_records<W>(writer: &mut W, header: &Header) -> io::Result<()>
where
    W: Write,
{
    let file_format = header.file_format();

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
        write_other(writer, file_format, key, collection)?;
    }

    Ok(())
}

/// Writes the records in the order they were added, followed by any record the header does not
/// track, grouped by category.
///
/// A key that no longer resolves to a record is skipped, which is what makes this tolerant of
/// headers mutated through, e.g., `Header::infos_mut`.
fn write_added_records<W>(writer: &mut W, header: &Header) -> io::Result<()>
where
    W: Write,
{
    let file_format = header.file_format();

    let mut infos = HashSet::new();
    let mut filters = HashSet::new();
    let mut formats = HashSet::new();
    let mut alternative_alleles = HashSet::new();
    let mut contigs = HashSet::new();
    let mut others: HashSet<&Other> = HashSet::new();

    for key in header.record_order() {
        match key {
            RecordKey::Info(id) => {
                if let Some(info) = header.infos().get(id.as_str()) {
                    write_info(writer, id, info)?;
                    infos.insert(id.as_str());
                }
            }
            RecordKey::Filter(id) => {
                if let Some(filter) = header.filters().get(id.as_str()) {
                    write_filter(writer, id, filter)?;
                    filters.insert(id.as_str());
                }
            }
            RecordKey::Format(id) => {
                if let Some(format) = header.formats().get(id.as_str()) {
                    write_format(writer, id, format)?;
                    formats.insert(id.as_str());
                }
            }
            RecordKey::AlternativeAllele(id) => {
                if let Some(alternative_allele) = header.alternative_alleles().get(id.as_str()) {
                    write_alternative_allele(writer, id, alternative_allele)?;
                    alternative_alleles.insert(id.as_str());
                }
            }
            RecordKey::Contig(id) => {
                if let Some(contig) = header.contigs().get(id.as_str()) {
                    write_contig(writer, id, contig)?;
                    contigs.insert(id.as_str());
                }
            }
            RecordKey::Other(key) => {
                if let Some(collection) = header.other_records().get(key) {
                    write_other(writer, file_format, key, collection)?;
                    others.insert(key);
                }
            }
        }
    }

    for (id, info) in header.infos() {
        if !infos.contains(id.as_str()) {
            write_info(writer, id, info)?;
        }
    }

    for (id, filter) in header.filters() {
        if !filters.contains(id.as_str()) {
            write_filter(writer, id, filter)?;
        }
    }

    for (id, format) in header.formats() {
        if !formats.contains(id.as_str()) {
            write_format(writer, id, format)?;
        }
    }

    for (id, alternative_allele) in header.alternative_alleles() {
        if !alternative_alleles.contains(id.as_str()) {
            write_alternative_allele(writer, id, alternative_allele)?;
        }
    }

    for (id, contig) in header.contigs() {
        if !contigs.contains(id.as_str()) {
            write_contig(writer, id, contig)?;
        }
    }

    for (key, collection) in header.other_records() {
        if !others.contains(key) {
            write_other(writer, file_format, key, collection)?;
        }
    }

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
        use crate::header::FileFormat;

        let mut buf = Vec::new();

        let header = Header::builder()
            .set_file_format(FileFormat::new(4, 5))
            .build();

        write_header(&mut buf, &header, RecordOrder::default())?;

        let expected = b"##fileformat=VCFv4.5
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(buf, expected);

        Ok(())
    }

    const INTERLEAVED_HEADER: &str = "##fileformat=VCFv4.5
##FILTER=<ID=q10,Description=\"Quality below 10\">
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

    #[test]
    fn test_write_header_with_record_order_as_added() -> io::Result<()> {
        let header: Header = INTERLEAVED_HEADER.parse().map_err(io::Error::other)?;

        let mut buf = Vec::new();
        write_header(&mut buf, &header, RecordOrder::AsAdded)?;

        assert_eq!(buf, INTERLEAVED_HEADER.as_bytes());

        Ok(())
    }

    #[test]
    fn test_write_header_with_record_order_grouped() -> io::Result<()> {
        let header: Header = INTERLEAVED_HEADER.parse().map_err(io::Error::other)?;

        let mut buf = Vec::new();
        write_header(&mut buf, &header, RecordOrder::Grouped)?;

        let expected = b"##fileformat=VCFv4.5
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples\">
##FILTER=<ID=q10,Description=\"Quality below 10\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_header_with_record_order_as_added_and_untracked_records() -> io::Result<()> {
        use crate::header::record::value::{Map, map::Filter};

        let mut header: Header = INTERLEAVED_HEADER.parse().map_err(io::Error::other)?;

        header
            .filters_mut()
            .insert(String::from("q20"), Map::<Filter>::new("Quality below 20"));

        let mut buf = Vec::new();
        write_header(&mut buf, &header, RecordOrder::AsAdded)?;

        let expected = b"##fileformat=VCFv4.5
##FILTER=<ID=q10,Description=\"Quality below 10\">
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples\">
##FILTER=<ID=q20,Description=\"Quality below 20\">
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
