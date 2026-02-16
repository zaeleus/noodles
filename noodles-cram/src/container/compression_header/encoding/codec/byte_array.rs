use std::{borrow::Cow, io};

use crate::{
    container::{
        block,
        compression_header::{
            Encoding,
            encoding::{
                Decode, Encode,
                codec::{Byte, Integer},
            },
        },
    },
    io::{
        BitReader, BitWriter, reader::container::slice::records::ExternalDataReaders,
        writer::container::slice::records::ExternalDataWriters,
    },
};

#[derive(Clone, Debug, Eq, PartialEq)]
#[allow(clippy::large_enum_variant)]
pub enum ByteArray {
    ByteArrayLength {
        len_encoding: Encoding<Integer>,
        value_encoding: Encoding<Byte>,
    },
    ByteArrayStop {
        stop_byte: u8,
        block_content_id: block::ContentId,
    },
}

impl<'de> Decode<'de> for ByteArray {
    type Value = Cow<'de, [u8]>;

    fn decode(
        &self,
        core_data_reader: &mut BitReader<'de>,
        external_data_readers: &mut ExternalDataReaders<'de>,
    ) -> std::io::Result<Self::Value> {
        match self {
            Self::ByteArrayLength {
                len_encoding,
                value_encoding,
            } => {
                let len = len_encoding
                    .decode(core_data_reader, external_data_readers)
                    .and_then(|n| {
                        usize::try_from(n)
                            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                    })?;

                value_encoding
                    .get()
                    .decode_take(core_data_reader, external_data_readers, len)
            }
            Self::ByteArrayStop {
                stop_byte,
                block_content_id,
            } => {
                let src = external_data_readers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                let Some(i) = src.iter().position(|&b| b == *stop_byte) else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "missing byte array stop byte",
                    ));
                };

                let (buf, rest) = src.split_at(i);
                *src = &rest[1..];

                Ok(Cow::Borrowed(buf))
            }
        }
    }
}

impl<'en> Encode<'en> for ByteArray {
    type Value = &'en [u8];

    fn encode(
        &self,
        core_data_writer: &mut BitWriter,
        external_data_writers: &mut ExternalDataWriters,
        value: Self::Value,
    ) -> io::Result<()> {
        match self {
            Self::ByteArrayLength {
                len_encoding,
                value_encoding,
            } => {
                let len = i64::try_from(value.len())
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                len_encoding.encode(core_data_writer, external_data_writers, len)?;

                for &v in value {
                    value_encoding.encode(core_data_writer, external_data_writers, v)?;
                }

                Ok(())
            }
            Self::ByteArrayStop {
                stop_byte,
                block_content_id,
            } => {
                let dst = external_data_writers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                dst.extend(value);
                dst.push(*stop_byte);

                Ok(())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode() -> io::Result<()> {
        fn t(
            external_data: &[u8],
            encoding: &Encoding<ByteArray>,
            expected: &[u8],
        ) -> io::Result<()> {
            let core_data = [];
            let mut core_data_reader = BitReader::new(&core_data[..]);

            let mut external_data_readers = ExternalDataReaders::new();
            external_data_readers.insert(1, external_data);

            let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;

            assert_eq!(expected, &*actual);

            Ok(())
        }

        let block_content_id = 1;
        let len_encoding = Encoding::new(Integer::External { block_content_id });
        let value_encoding = Encoding::new(Byte::External { block_content_id });
        t(
            &[0x04, 0x6e, 0x64, 0x6c, 0x73],
            &Encoding::new(ByteArray::ByteArrayLength {
                len_encoding,
                value_encoding,
            }),
            b"ndls",
        )?;

        t(
            &[0x6e, 0x64, 0x6c, 0x73, 0x00],
            &Encoding::new(ByteArray::ByteArrayStop {
                stop_byte: 0x00,
                block_content_id: 1,
            }),
            b"ndls",
        )?;

        assert!(matches!(
            t(
                &[0x6e, 0x64, 0x6c, 0x73],
                &Encoding::new(ByteArray::ByteArrayStop{ stop_byte: 0x00, block_content_id: 1 }),
                b""
            ),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_encode() -> io::Result<()> {
        fn t(
            encoding: &Encoding<ByteArray>,
            value: &[u8],
            expected_core_data: &[u8],
            expected_external_data: &[u8],
        ) -> io::Result<()> {
            let mut core_data_writer = BitWriter::default();

            let block_content_id = 1;
            let mut external_data_writers = [(block_content_id, Vec::new())].into_iter().collect();

            encoding.encode(&mut core_data_writer, &mut external_data_writers, value)?;

            let actual_core_data = core_data_writer.finish()?;
            assert_eq!(actual_core_data, expected_core_data);

            let actual_external_data = &external_data_writers[&block_content_id];
            assert_eq!(actual_external_data, expected_external_data);

            Ok(())
        }

        let block_content_id = 1;
        let len_encoding = Encoding::new(Integer::External { block_content_id });
        let value_encoding = Encoding::new(Byte::External { block_content_id });
        t(
            &Encoding::new(ByteArray::ByteArrayLength {
                len_encoding,
                value_encoding,
            }),
            b"ndls",
            &[],
            &[0x04, b'n', b'd', b'l', b's'],
        )?;

        t(
            &Encoding::new(ByteArray::ByteArrayStop {
                stop_byte: 0x00,
                block_content_id: 1,
            }),
            b"ndls",
            &[],
            &[b'n', b'd', b'l', b's', 0x00],
        )?;

        Ok(())
    }
}
