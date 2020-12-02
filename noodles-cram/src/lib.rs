mod bit_reader;
mod bit_writer;
pub mod container;
pub mod crai;
mod data_container;
pub mod file_definition;
mod huffman;
mod num;
mod rans;
pub mod reader;
pub mod record;
pub mod writer;

pub use self::{
    bit_reader::BitReader, bit_writer::BitWriter, container::Container,
    data_container::DataContainer, file_definition::FileDefinition, reader::Reader, record::Record,
    writer::Writer,
};

use std::{convert::TryFrom, fs::File, io, path::Path};

use noodles_bam as bam;

static MAGIC_NUMBER: &[u8] = b"CRAM";

pub fn index<P>(src: P) -> io::Result<crai::Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    reader.read_file_definition()?;
    reader.read_file_header()?;

    let mut index = Vec::new();
    let mut container_position = reader.position()?;

    loop {
        let container = reader.read_container()?;

        if container.is_eof() {
            break;
        }

        let container_header = container.header();
        let container_len = container_header.len();

        let landmarks = container_header.landmarks();
        let slice_count = landmarks.len();

        let data_container = DataContainer::try_from(container.clone())?;

        for (i, slice) in data_container.slices().iter().enumerate() {
            let slice_header = slice.header();

            let slice_reference_sequence_id = slice_header.reference_sequence_id();

            if slice_reference_sequence_id.is_many() {
                todo!("unhandled multi-reference slice");
            }

            let reference_sequence_id = if slice_reference_sequence_id.is_none() {
                None
            } else {
                bam::record::ReferenceSequenceId::try_from(i32::from(slice_reference_sequence_id))
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
            };

            let landmark = landmarks[i];

            let slice_length = if i < slice_count - 1 {
                landmarks[i + 1] - landmark
            } else {
                container_len - landmark
            };

            let record = crai::Record::new(
                reference_sequence_id,
                slice_header.alignment_start(),
                slice_header.alignment_span(),
                container_position,
                landmark as u64,
                slice_length as u64,
            );

            index.push(record);
        }

        container_position = reader.position()?;
    }

    Ok(index)
}
