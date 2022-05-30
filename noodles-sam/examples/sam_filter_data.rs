//! Filters records in a SAM that have a unique alignment.
//!
//! That is, there is a single alignment hit count (SAM record data tag `NM` = 1).

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_sam::{self as sam, record::data::field::Tag, AlignmentRecord};

fn is_unique_record(record: &sam::Record) -> io::Result<bool> {
    match record.data().get(Tag::AlignmentHitCount) {
        Some(field) => field.value().as_int().map(|hits| hits == 1).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected integer, got {:?}", field.value()),
            )
        }),
        None => Ok(false),
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");

    let mut reader = File::open(src).map(BufReader::new).map(sam::Reader::new)?;
    let header: sam::Header = reader.read_header()?.parse()?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = sam::Writer::new(handle);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;

        if is_unique_record(&record)? {
            writer.write_record(&header, &record)?;
        }
    }

    Ok(())
}
