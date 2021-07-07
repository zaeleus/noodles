//! Compute the depth over a set of positions.

use noodles_bam as bam;
use noodles_core::Region;
use std::{env, fs::File, io};
fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");
    let src_idx = env::args().nth(2).expect("missing src idx");

    let mut reader = File::open(src).map(bam::Reader::new)?;
    let _header = reader.read_header()?;

    let index = bam::bai::read(src_idx)?;
    let ref_seqs = reader.read_reference_sequences()?;

    let info = ref_seqs.get("chr1").unwrap();

    let pileups = bam::Pileup::new(reader.query(
        &ref_seqs,
        &index,
        &Region::mapped("chr1", 1, info.len()),
    )?)?;

    let mut count = 0;
    for pileup in pileups {
        println!("Depth: {:?}", pileup?.depth);
        count += 1;
    }
    println!("count: {}", count);

    Ok(())
}
