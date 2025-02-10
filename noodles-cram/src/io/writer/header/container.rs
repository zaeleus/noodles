use std::{
    io::{self, Write},
    mem,
};

use flate2::Compression;
use noodles_sam as sam;

use crate::{
    codecs::Encoder,
    container::{block::ContentType, Block},
    data_container::Header,
    io::writer::container::{write_block, write_header},
};

pub(super) fn write_container<W>(writer: &mut W, header: &sam::Header) -> io::Result<()>
where
    W: Write,
{
    const ENCODER: Encoder = Encoder::Gzip(Compression::new(6));

    validate_reference_sequences(header.reference_sequences())?;

    let data = serialize_header(header)?;

    let block = Block::builder()
        .set_content_type(ContentType::FileHeader)
        .compress_and_set_data(data, ENCODER)?
        .build();

    let header = build_header();
    write_header(writer, &header, block.len())?;

    write_block(writer, &block)?;

    Ok(())
}

fn validate_reference_sequences(
    reference_sequences: &sam::header::ReferenceSequences,
) -> io::Result<()> {
    use sam::header::record::value::map::reference_sequence::tag;

    for reference_sequence in reference_sequences.values() {
        if !reference_sequence
            .other_fields()
            .contains_key(&tag::MD5_CHECKSUM)
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "header reference sequence record is missing MD5 checksum",
            ));
        }
    }

    Ok(())
}

fn serialize_header(header: &sam::Header) -> io::Result<Vec<u8>> {
    const LENGTH_SIZE: usize = mem::size_of::<i32>();

    let mut buf = vec![0; LENGTH_SIZE];

    let mut writer = sam::io::Writer::new(&mut buf);
    writer.write_header(header)?;

    // SAFETY: `buf.len() >= LENGTH_SIZE`.
    let len = i32::try_from(buf.len() - LENGTH_SIZE)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    buf[..LENGTH_SIZE].copy_from_slice(&len.to_le_bytes());

    Ok(buf)
}

fn build_header() -> Header {
    Header::builder().set_block_count(1).build()
}

#[cfg(test)]
mod tests {
    use bytes::BufMut;

    use super::*;

    #[test]
    fn test_write_container() -> Result<(), Box<dyn std::error::Error>> {
        use byteorder::{LittleEndian, WriteBytesExt};
        use flate2::CrcWriter;
        use sam::header::record::value::{
            map::{self, header::Version},
            Map,
        };

        use crate::{codecs::gzip, io::writer::num::write_itf8};

        let header_header = Map::<map::Header>::new(Version::new(1, 6));
        let header = sam::Header::builder().set_header(header_header).build();

        let mut buf = Vec::new();
        write_container(&mut buf, &header)?;

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

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_validate_reference_sequences() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZeroUsize;

        use sam::header::record::value::{
            map::{reference_sequence::tag, ReferenceSequence},
            Map,
        };

        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        let header = sam::Header::builder()
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::builder()
                    .set_length(SQ0_LN)
                    .insert(
                        tag::MD5_CHECKSUM,
                        Vec::from("d7eba311421bbc9d3ada44709dd61534"),
                    )
                    .build()?,
            )
            .build();
        assert!(validate_reference_sequences(header.reference_sequences()).is_ok());

        let header = sam::Header::builder()
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(SQ0_LN))
            .build();
        assert!(validate_reference_sequences(header.reference_sequences()).is_err());

        Ok(())
    }
}
