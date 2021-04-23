use std::{
    convert::TryFrom,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::record::{Ids, QualityScore};

use super::value::{read_value, Value};

#[allow(dead_code, unused_variables)]
pub fn read_site<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    let chrom = reader.read_i32::<LittleEndian>()?;
    let pos = reader.read_i32::<LittleEndian>()?;

    let rlen = reader.read_i32::<LittleEndian>()?;

    let qual = reader.read_f32::<LittleEndian>().and_then(|value| {
        QualityScore::try_from(value).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let n_allele_info = reader.read_i32::<LittleEndian>()?;
    let allele_count = n_allele_info >> 16;
    let info_count = n_allele_info & 0xffff;

    let n_fmt_sample = reader.read_i32::<LittleEndian>()?;
    let format_count = n_fmt_sample >> 24;
    let sample_count = n_fmt_sample & 0xffffff;

    let id = match read_value(reader) {
        Ok(Value::String(Some(id))) => id
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        Ok(Value::String(None)) => Ids::default(),
        Ok(v) => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected string, got {:?}", v),
            ))
        }
        Err(e) => return Err(e),
    };

    let mut ref_alt = Vec::with_capacity(allele_count as usize);

    for _ in 0..allele_count {
        match read_value(reader) {
            Ok(Value::String(Some(s))) => ref_alt.push(s),
            Ok(Value::String(None)) => ref_alt.push(String::from(".")),
            Ok(v) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("expected string, got {:?}", v),
                ))
            }
            Err(e) => return Err(e),
        }
    }

    let filter_indices = match read_value(reader) {
        Ok(Value::Int8(None)) => Vec::new(),
        Ok(Value::Int8(Some(i))) => vec![i],
        Ok(Value::Int8Array(indicies)) => indicies,
        Ok(v) => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected i8, got {:?}", v),
            ))
        }
        Err(e) => return Err(e),
    };

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_site() -> io::Result<()> {
        // § Putting it all together (2021-01-13)
        //
        // Note that the data in the reference table mixes big and little endian.
        let data = [
            0x01, 0x00, 0x00, 0x00, // chrom = 1
            0x64, 0x00, 0x00, 0x00, // pos = 100 (base 0)
            0x01, 0x00, 0x00, 0x00, // rlen = 1
            0xcd, 0xcc, 0xf0, 0x41, // qual = 30.1
            0x04, 0x00, 0x02, 0x00, // n_allele_info (allele count, info count) = (2, 4)
            0x03, 0x00, 0x00, 0x05, // n_fmt_sample (format count, sample_count) = (5, 3)
            0x57, 0x72, 0x73, 0x31, 0x32, 0x33, // id = "rs123"
            0x17, 0x41, // ref = A
            0x17, 0x43, // alt = C
            0x11, 0x00, // filter = 0 (PASS)
        ];
        let mut reader = &data[..];

        read_site(&mut reader)
    }
}
