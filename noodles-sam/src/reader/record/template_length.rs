use std::io;

pub(crate) fn parse_template_length(src: &[u8]) -> io::Result<i32> {
    lexical_core::parse(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
