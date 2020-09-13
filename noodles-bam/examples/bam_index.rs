//! Builds and writes a BAM index from a BAM file.

use std::{env, fs::File, io};

use noodles_bam::{
    self as bam,
    bai::{self, index::reference_sequence::bin::Chunk},
};
use noodles_sam as sam;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bam::Reader::new)?;
    let header: sam::Header = reader.read_header()?.parse()?;
    reader.read_reference_sequences()?;

    let mut record = bam::Record::default();

    let mut builder = bai::Index::builder();
    let mut start_position = reader.virtual_position();

    loop {
        match reader.read_record(&mut record) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => return Err(e.into()),
        }

        let end_position = reader.virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        builder.add_record(&record, chunk);

        start_position = end_position;
    }

    let index = builder.build(header.reference_sequences().len());

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = bai::Writer::new(handle);

    writer.write_header()?;
    writer.write_index(&index)?;

    Ok(())
}
