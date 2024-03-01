mod keys;
mod sample;

use std::io::{self, Write};

use self::{keys::write_keys, sample::write_sample};
use crate::{
    variant::{record::Samples as _, record_buf::Samples},
    Header,
};

pub(super) fn write_samples<W>(writer: &mut W, header: &Header, samples: &Samples) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b"\t";

    write_keys(writer, samples.keys())?;

    for sample in samples.iter() {
        writer.write_all(DELIMITER)?;
        write_sample(writer, header, sample)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_samples() -> Result<(), Box<dyn std::error::Error>> {
        use crate::variant::record_buf::samples::{keys::key, sample::Value, Keys};

        fn t(
            buf: &mut Vec<u8>,
            header: &Header,
            genotypes: &Samples,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_samples(buf, header, genotypes)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let header = Header::default();
        let mut buf = Vec::new();

        let genotypes = Samples::new(
            Keys::try_from(vec![String::from(key::GENOTYPE)])?,
            vec![vec![Some(Value::from("0|0"))]],
        );
        t(&mut buf, &header, &genotypes, b"GT\t0|0")?;

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
        t(&mut buf, &header, &genotypes, b"GT:GQ\t0|0:13\t0/1:8")?;

        Ok(())
    }
}
