//! Filters records in a BAM that have a unique alignment.
//!
//! That is, there is a single alignment hit count (SAM record data tag `NM` = 1).

use std::{env, fs::File, io};

use noodles_bam as bam;
use noodles_sam::{self as sam, record::data::field::Tag};

fn is_unique_record(record: &bam::Record) -> io::Result<bool> {
    for result in record.data().fields() {
        let field = result?;

        if field.tag() == &Tag::AlignmentHitCount {
            let value = field.value();

            if let Some(hits) = value.as_int() {
                return Ok(hits == 1);
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("expected integer, got {:?}", value),
                ));
            }
        }
    }

    Ok(false)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");

    let mut reader = File::open(src).map(bam::Reader::new)?;
    let header: sam::Header = reader.read_header()?.parse()?;
    reader.read_reference_sequences()?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = bam::Writer::new(handle);

    writer.write_header(&header)?;
    writer.write_reference_sequences(header.reference_sequences())?;

    for result in reader.records() {
        let record = result?;

        if is_unique_record(&record)? {
            writer.write_record(&record)?;
        }
    }

    Ok(())
}
