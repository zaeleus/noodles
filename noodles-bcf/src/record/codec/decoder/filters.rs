use std::io;

use super::string_map::read_string_map_indices;

pub(crate) fn read_filter(reader: &mut &[u8]) -> io::Result<Vec<usize>> {
    read_string_map_indices(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
