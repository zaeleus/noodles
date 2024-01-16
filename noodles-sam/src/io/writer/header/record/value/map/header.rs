mod version;

use std::io::{self, Write};

use self::version::write_version_field;
use super::write_other_fields;
use crate::header::record::value::{map::Header, Map};

pub(crate) fn write_header<W>(writer: &mut W, header: &Map<Header>) -> io::Result<()>
where
    W: Write,
{
    write_version_field(writer, header.version())?;
    write_other_fields(writer, header.other_fields())?;
    Ok(())
}
