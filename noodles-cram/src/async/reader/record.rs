use std::{
    collections::HashMap,
    convert::TryFrom,
    io::{self, Read},
};

use noodles_sam as sam;
use tokio::io::{AsyncBufRead, AsyncRead};

use crate::{
    container::ReferenceSequenceId,
    data_container::{compression_header::Encoding, CompressionHeader},
    num::Itf8,
    r#async::reader::num::read_itf8,
    BitReader, Record,
};

pub struct Reader<'a, CDR, EDR>
where
    CDR: Read,
    EDR: AsyncBufRead + Unpin,
{
    compression_header: &'a CompressionHeader,
    core_data_reader: BitReader<CDR>,
    external_data_readers: HashMap<Itf8, EDR>,
    reference_sequence_id: ReferenceSequenceId,
    prev_alignment_start: Itf8,
}

impl<'a, CDR, EDR> Reader<'a, CDR, EDR>
where
    CDR: Read,
    EDR: AsyncBufRead + Unpin,
{
    pub fn new(
        compression_header: &'a CompressionHeader,
        core_data_reader: BitReader<CDR>,
        external_data_readers: HashMap<Itf8, EDR>,
        reference_sequence_id: ReferenceSequenceId,
        initial_alignment_start: Itf8,
    ) -> Self {
        Self {
            compression_header,
            core_data_reader,
            external_data_readers,
            reference_sequence_id,
            prev_alignment_start: initial_alignment_start,
        }
    }

    pub async fn read_record(&mut self) -> io::Result<Record> {
        let mut builder = Record::builder();

        let bam_bit_flags = self.read_bam_bit_flags().await?;
        builder = builder.set_bam_flags(bam_bit_flags);

        Ok(builder.build())
    }

    async fn read_bam_bit_flags(&mut self) -> io::Result<sam::record::Flags> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .bam_bit_flags_encoding();

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
        .and_then(|n| u16::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
        .map(sam::record::Flags::from)
    }
}

async fn decode_itf8<CDR, EDR>(
    encoding: &Encoding,
    _core_data_reader: &mut BitReader<CDR>,
    external_data_readers: &mut HashMap<Itf8, EDR>,
) -> io::Result<Itf8>
where
    CDR: Read,
    EDR: AsyncRead + Unpin,
{
    match encoding {
        Encoding::External(block_content_id) => {
            let reader = external_data_readers
                .get_mut(block_content_id)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "missing external block")
                })?;

            read_itf8(reader).await
        }
        _ => todo!("decode_itf8: {:?}", encoding),
    }
}
