pub(crate) mod data_series_encodings;
mod encoding;
mod preservation_map;
mod tag_encodings;

use std::io::{self, Write};

use self::{
    data_series_encodings::write_data_series_encodings,
    encoding::{write_byte_array_encoding, write_byte_encoding, write_integer_encoding},
    preservation_map::{build_preservation_map, write_preservation_map},
    tag_encodings::{build_tag_encodings, write_tag_encodings},
};
use crate::{
    container::{CompressionHeader, compression_header::DataSeriesEncodings},
    file_definition::Version,
    io::writer::{Options, Record},
};

pub fn write_compression_header<W>(
    writer: &mut W,
    compression_header: &CompressionHeader,
    version: Version,
) -> io::Result<()>
where
    W: Write,
{
    write_preservation_map(writer, compression_header.preservation_map(), version)?;
    write_data_series_encodings(writer, compression_header.data_series_encodings(), version)?;
    write_tag_encodings(writer, compression_header.tag_encodings(), version)?;
    Ok(())
}

pub(super) fn build_compression_header(options: &Options, records: &[Record]) -> CompressionHeader {
    CompressionHeader {
        preservation_map: build_preservation_map(options, records),
        data_series_encodings: DataSeriesEncodings::init(options.version),
        tag_encodings: build_tag_encodings(records, options.version),
    }
}
