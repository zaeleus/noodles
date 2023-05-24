use std::io;

use byteorder::WriteBytesExt;

use crate::{
    container::block,
    data_container::compression_header::{
        encoding::{
            codec::{Byte, Integer},
            Decode, Encode,
        },
        Encoding,
    },
    io::BitReader,
    reader::record::ExternalDataReaders,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ByteArray {
    // len_encoding, value_encoding
    ByteArrayLen(Encoding<Integer>, Encoding<Byte>),
    // stop_byte, block_content_id
    ByteArrayStop(u8, block::ContentId),
}

impl Decode for ByteArray {
    type Value = Vec<u8>;

    fn decode<R, S>(
        &self,
        core_data_reader: &mut BitReader<R>,
        external_data_readers: &mut ExternalDataReaders<S>,
    ) -> std::io::Result<Self::Value>
    where
        R: bytes::Buf,
        S: bytes::Buf,
    {
        match self {
            ByteArray::ByteArrayLen(len_encoding, value_encoding) => {
                let len = len_encoding.decode(core_data_reader, external_data_readers)?;

                let mut buf = vec![0; len as usize];

                for value in &mut buf {
                    *value = value_encoding.decode(core_data_reader, external_data_readers)?;
                }

                Ok(buf)
            }
            ByteArray::ByteArrayStop(stop_byte, block_content_id) => {
                let src = external_data_readers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                let len = match src.chunk().iter().position(|&b| b == *stop_byte) {
                    Some(i) => i,
                    None => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "missing byte array stop byte",
                        ))
                    }
                };

                let mut buf = vec![0; len];
                src.copy_to_slice(&mut buf);

                // Discard the stop byte.
                src.advance(1);

                Ok(buf)
            }
        }
    }
}

impl<'en> Encode<'en> for ByteArray {
    type Value = &'en [u8];

    fn encode<W, X>(
        &self,
        core_data_writer: &mut crate::io::BitWriter<W>,
        external_data_writers: &mut std::collections::HashMap<block::ContentId, X>,
        value: Self::Value,
    ) -> io::Result<()>
    where
        W: io::Write,
        X: io::Write,
    {
        match self {
            ByteArray::ByteArrayLen(len_encoding, value_encoding) => {
                let len = i32::try_from(value.len())
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                len_encoding.encode(core_data_writer, external_data_writers, len)?;

                for &v in value {
                    value_encoding.encode(core_data_writer, external_data_writers, v)?;
                }

                Ok(())
            }
            ByteArray::ByteArrayStop(stop_byte, block_content_id) => {
                let writer = external_data_writers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                writer.write_all(value)?;
                writer.write_u8(*stop_byte)?;

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
            external_data_readers.insert(block::ContentId::from(1), external_data);

            let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;

            assert_eq!(expected, actual);

            Ok(())
        }

        let len_encoding = Encoding::new(Integer::External(block::ContentId::from(1)));
        let value_encoding = Encoding::new(Byte::External(block::ContentId::from(1)));
        t(
            &[0x04, 0x6e, 0x64, 0x6c, 0x73],
            &Encoding::new(ByteArray::ByteArrayLen(len_encoding, value_encoding)),
            b"ndls",
        )?;

        t(
            &[0x6e, 0x64, 0x6c, 0x73, 0x00],
            &Encoding::new(ByteArray::ByteArrayStop(0x00, block::ContentId::from(1))),
            b"ndls",
        )?;

        assert!(matches!(
            t(
                &[0x6e, 0x64, 0x6c, 0x73],
                &Encoding::new(ByteArray::ByteArrayStop(0x00, block::ContentId::from(1))),
                b""
            ),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
