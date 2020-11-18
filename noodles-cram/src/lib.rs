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

        let container_reference_sequence_id = container_header.reference_sequence_id();

        if container_reference_sequence_id.is_many() {
            todo!("unhandled multi-reference slice");
        }

        let reference_sequence_id = if container_reference_sequence_id.is_none() {
            None
        } else {
            bam::record::ReferenceSequenceId::try_from(i32::from(container_reference_sequence_id))
                .map(Some)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
        };

        let alignment_start = container_header.start_position();
        let alignment_span = container_header.alignment_span();

        let landmarks = container_header.landmarks();

        if landmarks.len() != 1 {
            todo!("unhandled multi-slice container");
        }

        let landmark = landmarks.first().copied().expect("missing landmark");
        let slice_length = container_len - landmark;

        let record = crai::Record::new(
            reference_sequence_id,
            alignment_start,
            alignment_span,
            container_position,
            landmark as u64,
            slice_length as u64,
        );

        index.push(record);

        container_position = reader.position()?;
    }

    Ok(index)
}
