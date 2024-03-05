mod field;

use std::io::{self, Write};

use noodles_vcf::{self as vcf, header::StringMaps, variant::record::Info};

use self::field::write_field;

pub fn write_info<W, I>(
    writer: &mut W,
    header: &vcf::Header,
    string_maps: &StringMaps,
    info: I,
) -> io::Result<()>
where
    W: Write,
    I: Info,
{
    for result in info.iter(header) {
        let (key, value) = result?;
        write_field(writer, string_maps.strings(), key, value)?;
    }

    Ok(())
}
