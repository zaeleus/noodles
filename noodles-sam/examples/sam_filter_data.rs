//! Filters records in a SAM that have a unique alignment.
//!
//! That is, there is a single alignment hit count (SAM record data tag `NM` = 1).

use std::{env, io};

use noodles_sam::{
    self as sam,
    alignment::record::{data::field::tag, Data},
};

fn is_unique_record(record: &sam::Record) -> io::Result<bool> {
    match record.data().get(&tag::ALIGNMENT_HIT_COUNT).transpose()? {
        Some(value) => value.as_int().map(|hits| hits == 1).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected an integer, got {:?}", value.ty()),
            )
        }),
        None => Ok(false),
    }
}

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");

    let mut reader = sam::io::reader::Builder::default().build_from_path(src)?;
    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = sam::io::Writer::new(stdout);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;

        if is_unique_record(&record)? {
            writer.write_record(&header, &record)?;
        }
    }

    Ok(())
}
