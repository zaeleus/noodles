use std::io;

use noodles_vcf::record::Ids;

use super::read_value;
use crate::lazy::record::Value;

pub fn read_id(src: &mut &[u8]) -> io::Result<Ids> {
    match read_value(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))? {
        Some(Value::String(Some(id))) => id
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        Some(Value::String(None)) => Ok(Ids::default()),
        v => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid id: expected string, got {v:?}"),
        )),
    }
}
