mod field;

use std::io::{self, Write};

use noodles_vcf::{self as vcf, header::string_maps::StringStringMap};

use self::field::write_field;

pub fn write_info<W>(
    writer: &mut W,
    string_string_map: &StringStringMap,
    info: &vcf::variant::record_buf::Info,
) -> io::Result<()>
where
    W: Write,
{
    for (key, value) in info.as_ref() {
        write_field(writer, string_string_map, key, value.as_ref())?;
    }

    Ok(())
}
