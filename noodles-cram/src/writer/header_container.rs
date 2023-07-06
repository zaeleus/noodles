mod header;

use std::io::{self, Write};

use bytes::BufMut;
use flate2::Compression;
use noodles_sam as sam;

use self::header::write_header;
use super::container::write_block;
use crate::{
    codecs::Encoder,
    container::{block::ContentType, Block},
};

pub fn write_header_container<W>(writer: &mut W, header: &sam::Header) -> io::Result<()>
where
    W: Write,
{
    const ENCODER: Encoder = Encoder::Gzip(Compression::new(6));

    validate_reference_sequences(header.reference_sequences())?;

    let header_data = header.to_string().into_bytes();
    let header_data_len = i32::try_from(header_data.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let mut data = Vec::new();
    data.put_i32_le(header_data_len);
    data.extend_from_slice(&header_data);

    let block = Block::builder()
        .set_content_type(ContentType::FileHeader)
        .compress_and_set_data(data, ENCODER)?
        .build();

    write_header(writer, block.len())?;
    write_block(writer, &block)?;

    Ok(())
}

fn validate_reference_sequences(
    reference_sequences: &sam::header::ReferenceSequences,
) -> io::Result<()> {
    for reference_sequence in reference_sequences.values() {
        if reference_sequence.md5_checksum().is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "header reference sequence record is missing MD5 checksum",
            ));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_header_container() -> Result<(), Box<dyn std::error::Error>> {
        use byteorder::{LittleEndian, WriteBytesExt};
        use flate2::CrcWriter;
        use sam::header::record::value::{
            map::{self, header::Version},
            Map,
        };

        use crate::{codecs::gzip, writer::num::write_itf8};

        let header_header = Map::<map::Header>::new(Version::new(1, 6));
        let header = sam::Header::builder().set_header(header_header).build();

        let mut actual = Vec::new();
        write_header_container(&mut actual, &header)?;

        let header_data = b"@HD\tVN:1.6\n";
        let header_data_len = i32::try_from(header_data.len())?;

        let mut data = Vec::new();
        data.put_i32_le(header_data_len);
        data.extend(header_data);

        let compressed_data = gzip::encode(Compression::new(6), &data)?;

        let mut block_writer = CrcWriter::new(Vec::new());
        block_writer.write_all(&[
            0x01, // compression method = 1 (Gzip)
            0x00, // content type = 0 (FileHeader)
            0x00, // block content ID = 0
        ])?;
        write_itf8(&mut block_writer, i32::try_from(compressed_data.len())?)?; // compressed length
        write_itf8(&mut block_writer, i32::try_from(data.len())?)?; // uncompressed length
        block_writer.write_all(&compressed_data)?;

        let crc32 = block_writer.crc().sum();
        let mut expected_block = block_writer.into_inner();
        expected_block.put_u32_le(crc32); // crc32

        let mut header_writer = CrcWriter::new(Vec::new());
        header_writer.write_i32::<LittleEndian>(i32::try_from(expected_block.len())?)?; // length
        write_itf8(&mut header_writer, -1)?; // reference sequence ID = -1 (None)
        header_writer.write_all(&[
            0x00, // alignment start = 0
            0x00, // alignment span = 0
            0x00, // record count = 0
            0x00, // record counter = 0
            0x00, // base count = 0
            0x01, // block count = 1
            0x00, // landmarks.len = 0
        ])?;

        let crc32 = header_writer.crc().sum();
        let mut expected_header = header_writer.into_inner();
        expected_header.put_u32_le(crc32); // crc32

        let mut expected = expected_header;
        expected.extend(expected_block);

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_validate_reference_sequences() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZeroUsize;

        use sam::header::record::value::{
            map::{reference_sequence::Md5Checksum, ReferenceSequence},
            Map,
        };

        let header = sam::Header::builder()
            .add_reference_sequence(
                "sq0".parse()?,
                Map::<ReferenceSequence>::builder()
                    .set_length(NonZeroUsize::try_from(8)?)
                    .set_md5_checksum(Md5Checksum::from([
                        0xd7, 0xeb, 0xa3, 0x11, 0x42, 0x1b, 0xbc, 0x9d, 0x3a, 0xda, 0x44, 0x70,
                        0x9d, 0xd6, 0x15, 0x34,
                    ]))
                    .build()?,
            )
            .build();
        assert!(validate_reference_sequences(header.reference_sequences()).is_ok());

        let header = sam::Header::builder()
            .add_reference_sequence(
                "sq0".parse()?,
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(8)?),
            )
            .build();
        assert!(validate_reference_sequences(header.reference_sequences()).is_err());

        Ok(())
    }
}
