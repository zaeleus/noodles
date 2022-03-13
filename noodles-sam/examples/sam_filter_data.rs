//! Filters records in a SAM that have a unique alignment.
//!
//! That is, there is a single alignment hit count (SAM record data tag `NM` = 1).

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_sam::{
    self as sam,
    record::data::field::{Tag, Value},
    AlignmentRecord,
};

fn is_unique_record(record: &sam::Record) -> io::Result<bool> {
    use sam::record::data::field::value::Type;

    let value = record
        .data()
        .get(Tag::AlignmentHitCount)
        .map(|field| field.value());

    match value {
        Some(Value::Int32(n)) => Ok(*n == 1),
        Some(v) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("expected {:?}, got {:?}", Type::Int32, v),
        )),
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
            writer.write_record(&record)?;
        }
    }

    Ok(())
}
