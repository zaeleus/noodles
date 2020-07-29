//! Counts the number of records in a CRAM file.
//!
//! Note that this example counts each successfully read record and not by summing the number of
//! records field in each container or slice header.
//!
//! The result matches the output of `samtools view -c <src>`.

use std::{convert::TryFrom, env, fs::File, io};

use noodles_cram as cram;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(cram::Reader::new)?;
    reader.read_file_definition()?;
    reader.read_file_header()?;

    let mut n = 0;

    loop {
        let container = reader.read_container()?;

        if container.is_eof() {
            break;
        }

        let data_container = cram::DataContainer::try_from(container)?;

        for slice in data_container.slices() {
            let records = slice.records(data_container.compression_header())?;
            n += records.len();
        }
    }

    println!("{}", n);

    Ok(())
}
