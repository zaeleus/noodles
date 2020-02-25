use std::{
    collections::HashMap,
    convert::TryFrom,
    io::{self, BufRead, Read},
};

use byteorder::ReadBytesExt;

use crate::{
    data_series::DataSeries,
    encoding::{self, Encoding},
    huffman::CanonicalHuffmanDecoder,
    num::{read_itf8, Itf8},
    BitReader, CompressionHeader, Record, Tag,
};

pub struct Reader<'a, R, S>
where
    R: Read,
    S: BufRead,
{
    compression_header: &'a CompressionHeader,
    core_data_reader: R,
    external_data_readers: HashMap<Itf8, S>,
    reference_sequence_id: Itf8,
    prev_alignment_start: Itf8,
}

impl<'a, R, S> Reader<'a, R, S>
where
    R: Read,
    S: BufRead,
{
    pub fn new(
        compression_header: &'a CompressionHeader,
        core_data_reader: R,
        external_data_readers: HashMap<Itf8, S>,
        reference_sequence_id: Itf8,
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

    pub fn read_record(&mut self, record: &mut Record) -> io::Result<()> {
        record.bam_bit_flags = self.read_bam_bit_flags()?;
        record.cram_bit_flags = self.read_cram_bit_flags()?;

        self.read_positional_data(record)?;

        let preservation_map = self.compression_header.preservation_map();

        if preservation_map.read_names_included() {
            record.read_name = self.read_read_name()?;
        } else {
            todo!("generate name");
        }

        self.read_mate_data(record)?;
        self.read_tag_data(record)?;

        Ok(())
    }

    fn read_bam_bit_flags(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::BamBitFlags)
            .expect("missing BF");

        decode_itf8(
            &encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    fn read_cram_bit_flags(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::CramBitFlags)
            .expect("missing CF");

        decode_itf8(
            &encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    fn read_positional_data(&mut self, record: &mut Record) -> io::Result<()> {
        record.reference_id = if self.reference_sequence_id == -2 {
            self.read_reference_id()?
        } else {
            self.reference_sequence_id
        };

        record.read_length = self.read_read_length()?;

        let ap_data_series_delta = self
            .compression_header
            .preservation_map()
            .ap_data_series_delta();

        let alignment_start = self.read_alignment_start()?;

        record.alignment_start = if ap_data_series_delta {
            self.prev_alignment_start + alignment_start
        } else {
            alignment_start
        };

        record.read_group = self.read_read_group()?;

        Ok(())
    }

    fn read_reference_id(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReferenceId)
            .expect("missing RI");

        decode_itf8(
            &encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    fn read_read_length(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReadLengths)
            .expect("missing RL");

        decode_itf8(
            &encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    fn read_alignment_start(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::InSeqPositions)
            .expect("missing AP");

        decode_itf8(
            &encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    fn read_read_group(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReadGroups)
            .expect("missing RG");

        decode_itf8(
            &encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    pub fn read_read_name(&mut self) -> io::Result<Vec<u8>> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReadNames)
            .expect("missing RN");

        decode_byte_array(
            &encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
            None,
        )
    }

    pub fn read_mate_data(&mut self, record: &mut Record) -> io::Result<()> {
        let cram_bit_flags = record.cram_bit_flags();

        if cram_bit_flags.is_detached() {
            record.next_mate_bit_flags = self.read_next_mate_bit_flags()?;

            let preservation_map = self.compression_header.preservation_map();

            if !preservation_map.read_names_included() {
                record.read_name = self.read_read_name()?;
            }

            record.next_fragment_reference_sequence_id =
                self.read_next_fragment_reference_sequence_id()?;
            record.next_mate_alignment_start = self.read_next_mate_alignment_start()?;
            record.template_size = self.read_template_size()?;
        } else if cram_bit_flags.has_mate_downstream() {
            record.distance_to_next_fragment = self.read_distance_to_next_fragment()?;
        }

        Ok(())
    }

    pub fn read_next_mate_bit_flags(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::NextMateBitFlags)
            .expect("missing MF");

        decode_itf8(
            &encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    pub fn read_next_fragment_reference_sequence_id(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::NextFragmentReferenceSequenceId)
            .expect("missing NS");

        decode_itf8(
            &encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    pub fn read_next_mate_alignment_start(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::NextMateAlignmentStart)
            .expect("missing NP");

        decode_itf8(
            &encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    pub fn read_template_size(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::TemplateSize)
            .expect("missing TS");

        decode_itf8(
            &encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    pub fn read_distance_to_next_fragment(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::DistanceToNextFragment)
            .expect("missing NF");

        decode_itf8(
            &encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    fn read_tag_data(&mut self, record: &mut Record) -> io::Result<()> {
        let tag_line = self.read_tag_line()?;

        let preservation_map = self.compression_header.preservation_map();
        let tag_ids_dictionary = preservation_map.tag_ids_dictionary();
        let tag_keys = &tag_ids_dictionary[tag_line as usize];

        let tag_encoding_map = self.compression_header.tag_encoding_map();

        record.tags.clear();

        for key in tag_keys {
            let id = (key[0] as i32) << 16 | (key[1] as i32) << 8 | (key[2] as i32);
            let encoding = tag_encoding_map.get(&id).expect("missing tag encoding");

            let data = decode_byte_array(
                encoding,
                &mut self.core_data_reader,
                &mut self.external_data_readers,
                None,
            )?;

            let tag = Tag::new([key[0], key[1]], data);
            record.add_tag(tag);
        }

        Ok(())
    }

    pub fn read_tag_line(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::TagIds)
            .expect("missing TL");

        decode_itf8(
            &encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }
}

fn decode_itf8<R, S>(
    encoding: &Encoding,
    core_data_reader: &mut R,
    external_data_readers: &mut HashMap<Itf8, S>,
) -> io::Result<Itf8>
where
    R: Read,
    S: Read,
{
    match encoding.kind() {
        encoding::Kind::External => {
            let mut reader = encoding.args();
            let block_content_id = read_itf8(&mut reader)?;

            let reader = external_data_readers
                .get_mut(&block_content_id)
                .expect("could not find block");

            read_itf8(reader)
        }
        encoding::Kind::Huffman => {
            let mut reader = encoding.args();

            let alphabet_len = read_itf8(&mut reader)? as usize;
            let mut alphabet = Vec::with_capacity(alphabet_len);

            for _ in 0..alphabet_len {
                let symbol = read_itf8(&mut reader)?;
                alphabet.push(symbol);
            }

            let bit_lens_len = read_itf8(&mut reader)? as usize;
            let mut bit_lens = Vec::with_capacity(bit_lens_len);

            for _ in 0..bit_lens_len {
                let len = read_itf8(&mut reader)?;
                bit_lens.push(len);
            }

            let decoder = CanonicalHuffmanDecoder::new(&alphabet, &bit_lens);
            let mut reader = BitReader::new(core_data_reader);
            decoder.read(&mut reader)
        }
        _ => todo!(),
    }
}

fn decode_byte_array<R, S>(
    encoding: &Encoding,
    core_data_reader: &mut R,
    external_data_readers: &mut HashMap<Itf8, S>,
    buf: Option<Vec<u8>>,
) -> io::Result<Vec<u8>>
where
    R: Read,
    S: BufRead,
{
    match encoding.kind() {
        encoding::Kind::External => {
            let mut reader = encoding.args();
            let block_content_id = read_itf8(&mut reader)?;

            let reader = external_data_readers
                .get_mut(&block_content_id)
                .expect("could not find block");

            let mut buf = buf.expect("missing buf");
            reader.read_exact(&mut buf)?;

            Ok(buf)
        }
        encoding::Kind::ByteArrayLen => {
            let mut reader = encoding.args();

            let len_encoding = read_encoding(&mut reader)?;
            let value_encoding = read_encoding(&mut reader)?;

            let len = decode_itf8(&len_encoding, core_data_reader, external_data_readers)?;

            let buf = vec![0; len as usize];
            let value = decode_byte_array(
                &value_encoding,
                core_data_reader,
                external_data_readers,
                Some(buf),
            )?;

            Ok(value)
        }
        encoding::Kind::ByteArrayStop => {
            let mut reader = encoding.args();
            let stop_byte = reader.read_u8()?;
            let block_content_id = read_itf8(&mut reader)?;

            let reader = external_data_readers
                .get_mut(&block_content_id)
                .expect("could not find block");

            let mut buf = Vec::new();
            reader.read_until(stop_byte, &mut buf)?;

            Ok(buf)
        }
        _ => todo!(),
    }
}

fn read_encoding<R>(reader: &mut R) -> io::Result<Encoding>
where
    R: Read,
{
    let kind = read_itf8(reader)
        .map(|codec_id| encoding::Kind::try_from(codec_id).expect("invalid codec id"))?;

    let args_len = read_itf8(reader)?;
    let mut args_buf = vec![0; args_len as usize];
    reader.read_exact(&mut args_buf)?;

    Ok(Encoding::new(kind, args_buf))
}
