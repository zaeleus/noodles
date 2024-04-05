mod field;

use std::io::{self, Write};

use self::field::write_field;
use super::MISSING;
use crate::{variant::record::Info, Header};

pub(super) fn write_info<W, I>(writer: &mut W, header: &Header, info: I) -> io::Result<()>
where
    W: Write,
    I: Info,
{
    const DELIMITER: &[u8] = b";";

    if info.is_empty() {
        writer.write_all(MISSING)?;
    } else {
        for (i, result) in info.iter(header).enumerate() {
            let (key, value) = result?;

            if i > 0 {
                writer.write_all(DELIMITER)?;
            }

            write_field(writer, key, value.as_ref())?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_info() -> io::Result<()> {
        use crate::variant::{
            record::info::field::key,
            record_buf::{info::field::Value as ValueBuf, Info as InfoBuf},
        };

        fn t(
            buf: &mut Vec<u8>,
            header: &Header,
            info: &InfoBuf,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_info(buf, header, info)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let header = Header::default();
        let mut buf = Vec::new();

        let info = InfoBuf::default();
        t(&mut buf, &header, &info, b".")?;

        let info = [(
            String::from(key::SAMPLES_WITH_DATA_COUNT),
            Some(ValueBuf::from(2)),
        )]
        .into_iter()
        .collect();
        t(&mut buf, &header, &info, b"NS=2")?;

        let info = [
            (
                String::from(key::SAMPLES_WITH_DATA_COUNT),
                Some(ValueBuf::from(2)),
            ),
            (String::from(key::IS_IN_DB_SNP), Some(ValueBuf::Flag)),
        ]
        .into_iter()
        .collect();

        t(&mut buf, &header, &info, b"NS=2;DB")?;

        Ok(())
    }
}
