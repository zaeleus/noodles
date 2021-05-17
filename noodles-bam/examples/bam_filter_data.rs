//! Filters records in a BAM that have a unique alignment.
//!
//! That is, there is a single alignment hit count (SAM record data tag `NM` = 1).

use std::{env, fs::File, io};

use noodles_bam::{self as bam, record::data::field::Value};
use noodles_sam::{self as sam, record::data::field::Tag};

fn is_unique_record(record: &bam::Record) -> io::Result<bool> {
    for result in record.data().fields() {
        let field = result?;

        if field.tag() == &Tag::AlignmentHitCount {
            let hits = match field.value() {
                Value::Int8(n) => i64::from(*n),
                Value::UInt8(n) => i64::from(*n),
                Value::Int16(n) => i64::from(*n),
                Value::UInt16(n) => i64::from(*n),
                Value::Int32(n) => i64::from(*n),
                Value::UInt32(n) => i64::from(*n),
                v => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "expected {{Int8, UInt8, Int16, UInt16, Int32, UInt32}}, got {:?}",
                            v
                        ),
                    ))
                }
            };

            return Ok(hits == 1);
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
