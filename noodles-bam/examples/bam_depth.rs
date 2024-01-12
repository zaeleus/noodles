//! Calculates the read depth of each position in a region.
//!
//! The results match the output of `samtools depth -r <region> <src>`.

use std::{
    env,
    io::{self, BufWriter, Write},
};

use noodles_bam as bam;
use noodles_core::Region;
use noodles_sam as sam;
use noodles_sam::alignment::iter::Depth;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");
    let region: Region = args.next().expect("missing region").parse()?;

    let mut reader = bam::io::indexed_reader::Builder::default().build_from_path(src)?;
    let header = reader.read_header()?;

    let query = reader
        .query(&header, &region)?
        .map(|result| result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>));

    let pileup = Depth::new(&header, query);

    let stdout = io::stdout().lock();
    let mut writer = BufWriter::new(stdout);

    let reference_sequence_name = region.name();

    for result in pileup {
        let (position, depth) = result?;

        if !region.interval().contains(position) {
            continue;
        }

        writeln!(writer, "{reference_sequence_name}\t{position}\t{depth}")?;
    }

    Ok(())
}
