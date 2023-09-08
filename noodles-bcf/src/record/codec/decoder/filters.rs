use std::io;

use crate::lazy::record::Filters;

pub(crate) fn read_filter(reader: &mut &[u8], filters: &mut Filters) -> io::Result<()> {
    use super::string_map::read_string_map_indices;

    let filter = filters.as_mut();
    filter.clear();

    let indices = read_string_map_indices(reader)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    filter.extend_from_slice(&indices);

    Ok(())
}
