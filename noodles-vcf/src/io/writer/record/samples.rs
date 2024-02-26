mod keys;
mod sample;

use std::io::{self, Write};

use self::{keys::write_keys, sample::write_sample};
use crate::variant::record_buf::Samples;

pub(super) fn write_samples<W>(writer: &mut W, genotypes: &Samples) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b"\t";

    write_keys(writer, genotypes.keys())?;

    for sample in genotypes.values() {
        writer.write_all(DELIMITER)?;
        write_sample(writer, &sample)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_samples() -> Result<(), Box<dyn std::error::Error>> {
        use crate::variant::record_buf::samples::{keys::key, sample::Value, Keys};

        fn t(buf: &mut Vec<u8>, genotypes: &Samples, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_samples(buf, genotypes)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let genotypes = Samples::new(
            Keys::try_from(vec![String::from(key::GENOTYPE)])?,
            vec![vec![Some(Value::from("0|0"))]],
        );
        t(&mut buf, &genotypes, b"GT\t0|0")?;

        let genotypes = Samples::new(
            Keys::try_from(vec![
                String::from(key::GENOTYPE),
                String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
            ])?,
            vec![
                vec![Some(Value::from("0|0")), Some(Value::from(13))],
                vec![Some(Value::from("0/1")), Some(Value::from(8))],
            ],
        );
        t(&mut buf, &genotypes, b"GT:GQ\t0|0:13\t0/1:8")?;

        Ok(())
    }
}
