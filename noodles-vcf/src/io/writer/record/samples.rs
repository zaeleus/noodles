mod keys;
mod sample;

use std::{
    error, fmt,
    io::{self, Write},
};

use self::{keys::write_keys, sample::write_sample};
use crate::{Header, variant::record::Samples};

/// An error returned when samples fail to write.
#[derive(Debug)]
pub enum WriteError {
    /// I/O error.
    Io(io::Error),
    /// The keys are invalid.
    InvalidKeys(keys::WriteError),
    /// A sample is invalid.
    InvalidSample(sample::WriteError),
    /// The sample count does not match the number of samples in the header.
    SampleCountMismatch {
        /// The actual sample count.
        actual: usize,
        /// The expected sample count.
        expected: usize,
    },
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidKeys(e) => Some(e),
            Self::InvalidSample(e) => Some(e),
            Self::SampleCountMismatch { .. } => None,
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidKeys(_) => write!(f, "invalid keys"),
            Self::InvalidSample(_) => write!(f, "invalid sample"),
            Self::SampleCountMismatch { actual, expected } => write!(
                f,
                "sample count mismatch: expected {expected}, got {actual}"
            ),
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

    let mut actual_sample_count = 0;

    for sample in samples.iter() {
        write_separator(writer)?;
        write_sample(writer, header, sample).map_err(WriteError::InvalidSample)?;
        actual_sample_count += 1;
    }

    let expected_sample_count = header.sample_names().len();

    if actual_sample_count != expected_sample_count {
        return Err(WriteError::SampleCountMismatch {
            actual: actual_sample_count,
            expected: expected_sample_count,
        });
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

        let mut buf = Vec::new();

        let header = Header::builder().add_sample_name("sample0").build();
        let genotypes = SamplesBuf::new(
            [String::from(key::GENOTYPE)].into_iter().collect(),
            vec![vec![Some(Value::from("0|0"))]],
        );
        t(&mut buf, &header, &genotypes, b"GT\t0|0")?;

        let header = Header::builder()
            .add_sample_name("sample0")
            .add_sample_name("sample1")
            .build();
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

        let header = Header::builder()
            .add_sample_name("sample0")
            .add_sample_name("sample1")
            .build();
        let genotypes = SamplesBuf::new(
            [String::from(key::GENOTYPE)].into_iter().collect(),
            vec![vec![Some(Value::from("0|0"))], vec![]],
        );
        t(&mut buf, &header, &genotypes, b"GT\t0|0\t.")?;

        buf.clear();
        let header = Header::default();
        let genotypes = SamplesBuf::new(
            [String::from(key::GENOTYPE)].into_iter().collect(),
            vec![vec![Some(Value::from("0|0"))]],
        );
        assert!(matches!(
            write_samples(&mut buf, &header, &genotypes),
            Err(WriteError::SampleCountMismatch {
                actual: 1,
                expected: 0
            })
        ));

        Ok(())
    }
}
