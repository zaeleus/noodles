//! Prints a CRAM file in the SAM format.
//!
//! The result matches the output of `samtools view <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
    path::Path,
};

use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

fn read_reference_assembly<P>(src: P) -> io::Result<Vec<fasta::Record>>
where
    P: AsRef<Path>,
{
    File::open(src)
        .map(BufReader::new)
        .map(fasta::Reader::new)?
        .records()
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let fasta_src = args.next().expect("missing fasta_src");
    let src = args.next().expect("missing src");

    let reference_assembly = read_reference_assembly(fasta_src)?;

    let mut reader = File::open(src).map(cram::Reader::new)?;
    reader.read_file_definition()?;

    let header: sam::Header = reader.read_file_header()?.parse()?;

    while let Some(container) = reader.read_data_container()? {
        for slice in container.slices() {
            for record in slice.records(container.compression_header())? {
                let sam_record = record.try_into_sam_record(
                    &reference_assembly,
                    header.reference_sequences(),
                    container.compression_header(),
                )?;

                println!("{}", sam_record);
            }
        }
    }

    Ok(())
}
