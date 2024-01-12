//! Concatenates BAM files.
//!
//! The result is similar to the output of `samtools cat --no-PG <srcs...>`.

use std::{env, io};

use noodles_bam as bam;

fn main() -> io::Result<()> {
    let srcs: Vec<_> = env::args().skip(1).collect();

    let first_src = srcs.first().expect("missing srcs[0]");
    let header = bam::io::reader::Builder
        .build_from_path(first_src)
        .and_then(|mut reader| reader.read_header())?;

    let stdout = io::stdout().lock();
    let mut writer = bam::io::Writer::new(stdout);

    writer.write_header(&header)?;

    for src in srcs {
        let mut reader = bam::io::reader::Builder.build_from_path(src)?;
        reader.read_header()?;

        io::copy(reader.get_mut(), writer.get_mut())?;
    }

    Ok(())
}
