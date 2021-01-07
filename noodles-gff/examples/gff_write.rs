//! Creates a new GFF file.
//!
//! This writes a GFF version directive and one (sparse) record to stdout.

use std::io;

use noodles_gff as gff;

fn main() -> io::Result<()> {
    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = gff::Writer::new(handle);

    let version = gff::Directive::GffVersion(Default::default());
    writer.write_directive(&version)?;

    let record = gff::Record::default();
    writer.write_record(&record)?;

    Ok(())
}
