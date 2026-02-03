mod keys;
mod sample;

use std::{
    error, fmt,
    io::{self, Write},
};

use self::{keys::write_keys, sample::write_sample};
use crate::{Header, variant::record::Samples};

/// An error returns when samples fail to write.
#[derive(Debug)]
pub enum WriteError {
    // I/O error.
    Io(io::Error),
    /// The keys are invalid.
    InvalidKeys(keys::WriteError),
    /// A sample is invalid.
    InvalidSample(sample::WriteError),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidKeys(e) => Some(e),
            Self::InvalidSample(e) => Some(e),
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidKeys(_) => write!(f, "invalid keys"),
            Self::InvalidSample(_) => write!(f, "invalid sample"),
        }
    }
}

pub(super) fn write_samples<W, S>(
    writer: &mut W,
    header: &Header,
    samples: S,
) -> Result<(), WriteError>
where
    W: Write,
    S: Samples,
{
    write_keys(writer, samples.column_names(header)).map_err(WriteError::InvalidKeys)?;

    for sample in samples.iter() {
        write_separator(writer)?;
        write_sample(writer, header, sample).map_err(WriteError::InvalidSample)?;
    }

    Ok(())
}

fn write_separator<W>(writer: &mut W) -> Result<(), WriteError>
where
    W: Write,
{
    const SEPARATOR: &[u8] = b"\t";
    writer.write_all(SEPARATOR).map_err(WriteError::Io)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record_buf::Samples as SamplesBuf;

    #[test]
    fn test_write_samples() -> Result<(), WriteError> {
        use crate::variant::{record::samples::keys::key, record_buf::samples::sample::Value};

        fn t(
            buf: &mut Vec<u8>,
            header: &Header,
            samples: &SamplesBuf,
            expected: &[u8],
        ) -> Result<(), WriteError> {
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
