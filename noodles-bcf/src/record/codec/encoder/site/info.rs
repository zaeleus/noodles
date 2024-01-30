mod field;

use std::io::{self, Write};

use noodles_vcf as vcf;

use self::field::write_field;
use crate::header::string_maps::StringStringMap;

pub fn write_info<W>(
    writer: &mut W,
    string_string_map: &StringStringMap,
    info: &vcf::record::Info,
) -> io::Result<()>
where
    W: Write,
{
    for (key, value) in info.as_ref() {
        write_field(writer, string_string_map, key, value.as_ref())?;
    }

    Ok(())
}
