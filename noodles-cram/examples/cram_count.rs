//! Counts the number of records in a CRAM file.
//!
//! Note that this example counts each successfully read record and not by summing the number of
//! records field in each container or slice header.
//!
//! The result matches the output of `samtools view --count <src>`.

use std::{env, io};

use noodles_cram::{self as cram, io::reader::Container};
use noodles_fasta as fasta;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = cram::io::reader::Builder::default().build_from_path(src)?;
    let header = reader.read_header()?;

    let reference_sequence_repository = fasta::Repository::default();

    let mut container = Container::default();
    let mut n = 0;

    while reader.read_container(&mut container)? != 0 {
        let compression_header = container.compression_header()?;

        for result in container.slices() {
            let slice = result?;

            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

            let records = slice.records(
                reference_sequence_repository.clone(),
                &header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            )?;

            n += records.len();
        }
    }

    println!("{n}");

    Ok(())
}
