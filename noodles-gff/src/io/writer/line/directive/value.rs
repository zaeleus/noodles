mod genome_build;
mod gff_version;
mod sequence_region;

use std::io::{self, Write};

use bstr::BStr;

pub(super) use self::{
    genome_build::write_genome_build, gff_version::write_gff_version,
    sequence_region::write_sequence_region,
};

pub(super) fn write_value<W>(writer: &mut W, value: &BStr) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(value)
}
