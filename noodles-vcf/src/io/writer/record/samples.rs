mod keys;
mod sample;

use std::io::{self, Write};

use self::{keys::write_keys, sample::write_sample};
use crate::{Header, variant::record::Samples};

pub(super) fn write_samples<W, S>(writer: &mut W, header: &Header, samples: S) -> io::Result<()>
where
    W: Write,
    S: Samples,
{
    const DELIMITER: &[u8] = b"\t";

    write_keys(writer, samples.column_names(header))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    for sample in samples.iter() {
        writer.write_all(DELIMITER)?;
        write_sample(writer, header, sample)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record_buf::Samples as SamplesBuf;

    #[test]
    fn test_write_samples() -> Result<(), Box<dyn std::error::Error>> {
        use crate::variant::{record::samples::keys::key, record_buf::samples::sample::Value};

        fn t(
            buf: &mut Vec<u8>,
            header: &Header,
            samples: &SamplesBuf,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_samples(buf, header, samples)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let header = Header::default();
        let mut buf = Vec::new();

        let genotypes = SamplesBuf::new(
            [String::from(key::GENOTYPE)].into_iter().collect(),
            vec![vec![Some(Value::from("0|0"))]],
        );
        t(&mut buf, &header, &genotypes, b"GT\t0|0")?;

        let genotypes = SamplesBuf::new(
            [
                String::from(key::GENOTYPE),
                String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
            ]
            .into_iter()
            .collect(),
            vec![
                vec![Some(Value::from("0|0")), Some(Value::from(13))],
                vec![Some(Value::from("0/1")), Some(Value::from(8))],
            ],
        );
        t(&mut buf, &header, &genotypes, b"GT:GQ\t0|0:13\t0/1:8")?;

        Ok(())
    }
}
