//! Creates a new GFF file.
//!
//! This writes a GFF version directive and one (sparse) record to stdout.

use std::io;

use noodles_gff::{self as gff, Line};

fn main() -> io::Result<()> {
    let stdout = io::stdout().lock();
    let mut writer = gff::io::Writer::new(stdout);

    let version = gff::Directive::GffVersion(Default::default());
    writer.write_directive(&version)?;

    let comment = Line::Comment(String::from("format: gff3"));
    writer.write_line(&comment)?;

    let record = gff::Record::default();
    writer.write_record(&record)?;

    Ok(())
}
