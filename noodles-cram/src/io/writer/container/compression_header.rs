pub(crate) mod data_series_encodings;
mod encoding;
mod preservation_map;
mod tag_encodings;

use std::io::{self, Write};

use self::{
    data_series_encodings::write_data_series_encodings,
    encoding::{
        write_encoding_for_byte_array_codec, write_encoding_for_byte_codec,
        write_encoding_for_integer_codec,
    },
    preservation_map::{build_preservation_map, write_preservation_map},
    tag_encodings::{build_tag_encodings, write_tag_encodings},
};
use crate::{
    container::{compression_header::DataSeriesEncodings, CompressionHeader},
    io::writer::{Options, Record},
};

pub fn write_compression_header<W>(
    writer: &mut W,
    compression_header: &CompressionHeader,
) -> io::Result<()>
where
    W: Write,
{
    write_preservation_map(writer, compression_header.preservation_map())?;
    write_data_series_encodings(writer, compression_header.data_series_encodings())?;
    write_tag_encodings(writer, compression_header.tag_encodings())?;
    Ok(())
}

pub(super) fn build_compression_header(options: &Options, records: &[Record]) -> CompressionHeader {
    CompressionHeader {
        preservation_map: build_preservation_map(options, records),
        data_series_encodings: DataSeriesEncodings::default(),
        tag_encodings: build_tag_encodings(records),
    }
}
