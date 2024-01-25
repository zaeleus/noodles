mod chromosome;
mod filters;
mod genotypes;
mod ids;
mod info;
mod quality_score;

use std::io::{self, Write};

use self::{
    chromosome::write_chromosome, filters::write_filters, genotypes::write_genotypes,
    ids::write_ids, info::write_info, quality_score::write_quality_score,
};
use crate::Record;

const MISSING: &[u8] = b".";

pub(super) fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b"\t";

    write_chromosome(writer, record.chromosome())?;

    writer.write_all(DELIMITER)?;
    write!(writer, "{}", record.position())?;

    writer.write_all(DELIMITER)?;
    write_ids(writer, record.ids())?;

    writer.write_all(DELIMITER)?;
    write!(writer, "{}", record.reference_bases())?;

    writer.write_all(DELIMITER)?;

    if record.alternate_bases().is_empty() {
        writer.write_all(MISSING)?;
    } else {
        write!(writer, "{}", record.alternate_bases())?;
    }

    writer.write_all(DELIMITER)?;
    write_quality_score(writer, record.quality_score())?;

    writer.write_all(DELIMITER)?;
    write_filters(writer, record.filters())?;

    writer.write_all(DELIMITER)?;
    write_info(writer, record.info())?;

    if !record.genotypes().is_empty() {
        writer.write_all(DELIMITER)?;
        write_genotypes(writer, record.genotypes())?;
    }

    writer.write_all(b"\n")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::Position;

    #[test]
    fn test_write_record() -> Result<(), Box<dyn std::error::Error>> {
        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::from(1))
            .set_reference_bases("A".parse()?)
            .build()?;

        let mut buf = Vec::new();
        write_record(&mut buf, &record)?;
        assert_eq!(buf, b"sq0\t1\t.\tA\t.\t.\t.\t.\n");

        Ok(())
    }
}
