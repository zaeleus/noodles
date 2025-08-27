//! Filters records in a BAM that have a unique alignment.
//!
//! That is, there is a single alignment hit count (SAM record data tag `NM` = 1).

use std::{env, fs::File, io};

use noodles_bam as bam;

fn is_unique_record(record: &bam::Record) -> io::Result<bool> {
    use noodles_sam::alignment::record::data::field::Tag;

    match record.data().get(&Tag::ALIGNMENT_HIT_COUNT).transpose()? {
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

    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = bam::io::Writer::new(stdout);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;

        if is_unique_record(&record)? {
            writer.write_record(&header, &record)?;
        }
    }

    Ok(())
}
