//! Rewrites a variant format to another variant format.
//!
//! The output format is determined from the extension of the destination.

use std::{env, io};

use noodles_util::variant;

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let dst = args.next().expect("missing dst");

    let mut reader = variant::io::reader::Builder::default().build_from_path(src)?;

    let header = reader.read_header()?;

    let mut writer = variant::io::writer::Builder::default().build_from_path(dst)?;

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    Ok(())
}
