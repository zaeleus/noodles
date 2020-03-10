//! Counts the number of records in a CRAM file.
//!
//! Note that this example counts each successfully read record and not by summing the number of
//! records field in each container or slice header.
//!
//! The result matches the output of `samtools view -c <src>`.

use std::{env, fs::File};

use noodles_cram as cram;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let file = File::open(src)?;
    let mut reader = cram::Reader::new(file);

    reader.read_file_definition()?;
    reader.read_file_header()?;

    let mut container = cram::Container::default();
    let mut n = 0;

    loop {
        container.clear();
        reader.read_container(&mut container)?;

        if container.is_eof() {
            break;
        }

        for slice in container.slices() {
            let mut record_reader = slice.records(container.compression_header());

            for _ in 0..slice.header().n_records() {
                let mut record = cram::Record::default();
                record_reader.read_record(&mut record)?;
                n += 1;
            }
        }
    }

    println!("{}", n);

    Ok(())
}
