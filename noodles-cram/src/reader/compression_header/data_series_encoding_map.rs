use std::{
    convert::TryFrom,
    io::{self, Read},
};

use crate::{
    container::compression_header::{
        data_series_encoding_map::DataSeries, encoding, DataSeriesEncodingMap, Encoding,
    },
    num::read_itf8,
};

pub fn read_data_series_encoding_map<R>(reader: &mut R) -> io::Result<DataSeriesEncodingMap>
where
    R: Read,
{
    let data_len = read_itf8(reader)?;
    let mut buf = vec![0; data_len as usize];
    reader.read_exact(&mut buf)?;

    let mut buf_reader = &buf[..];
    let map_len = read_itf8(&mut buf_reader)?;

    let mut encodings = DataSeriesEncodingMap::default();
    let mut key_buf = [0; 2];

    for _ in 0..map_len {
        buf_reader.read_exact(&mut key_buf)?;

        let key = DataSeries::try_from(&key_buf[..])
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let kind = read_itf8(&mut buf_reader).and_then(|codec_id| {
            encoding::Kind::try_from(codec_id)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let args_len = read_itf8(&mut buf_reader)?;
        let mut args_buf = vec![0; args_len as usize];
        buf_reader.read_exact(&mut args_buf)?;

        let encoding = Encoding::new(kind, args_buf);

        encodings.insert(key, encoding);
    }

    Ok(encodings)
}
