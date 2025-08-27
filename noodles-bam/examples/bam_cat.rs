//! Concatenates BAM files.
//!
//! The result is similar to the output of `samtools cat --no-PG <srcs...>`.

use std::{env, fs::File, io};

use noodles_bam as bam;

fn main() -> io::Result<()> {
    let srcs: Vec<_> = env::args().skip(1).collect();

    let first_src = srcs.first().expect("missing srcs[0]");
    let header = File::open(first_src)
        .map(bam::io::Reader::new)
        .and_then(|mut reader| reader.read_header())?;

    let stdout = io::stdout().lock();
    let mut writer = bam::io::Writer::new(stdout);

    writer.write_header(&header)?;

    for src in srcs {
        let mut reader = File::open(src).map(bam::io::Reader::new)?;
        reader.read_header()?;

        io::copy(reader.get_mut(), writer.get_mut())?;
    }

    Ok(())
}
