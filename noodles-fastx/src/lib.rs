pub mod error;
#[warn(missing_docs)]
/* mod declaration */
pub mod reader;
pub mod record;
pub mod records;
pub mod writer;

/* pub use declaration*/
pub use reader::Reader;
pub use record::Record;
pub use records::Records;
pub use writer::Writer;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_write() -> std::io::Result<()> {
        let data = b"\
@noodles:1/1
AGCT
+
abcd
@noodles: 2/1
TCGA
+noodles:2/1
dcba
>fasta sequence in middle of my fastq
GATC
ATGA
>aue
AACAC
@2
AGCT
+strange description
!!!!!!!!!
";

        let mut reader = Reader::new(&data[..])?;
        let mut writer = Writer::new(Vec::new());

        let record = reader.next_record().unwrap()?;
        assert!(record.is_fastq());
        writer.write_record(&record)?;

        let record = reader.next_record().unwrap()?;
        assert!(record.is_fastq());
        writer.write_record(&record)?;

        let record = reader.next_record().unwrap()?;
        assert!(!record.is_fastq());
        writer.write_record(&record)?;
        let record = reader.next_record().unwrap()?;
        assert!(!record.is_fastq());
        writer.write_record(&record)?;

        let mut record = reader.next_record().unwrap()?;
        assert!(record.is_fastq());
        record.fastq2fasta();
        assert!(!record.is_fastq());
        writer.write_record(&record)?;

        assert!(reader.next_record().is_none());

        let expected_result = b"\
@noodles:1/1
AGCT
+
abcd
@noodles: 2/1
TCGA
+noodles:2/1
dcba
>fasta sequence in middle of my fastq
GATCATGA
>aue
AACAC
>2
AGCT
";

        assert_eq!(
            String::from_utf8(expected_result.to_vec()).unwrap(),
            String::from_utf8(writer.get_ref().to_vec()).unwrap()
        );

        Ok(())
    }
}
