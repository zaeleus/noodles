//! Convert a Fastq file in Fasta

use std::{
    env,
    fs::File,
    io::stdout,
    io::{self, BufReader},
};

use noodles_fastx as fastx;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(|x| BufReader::with_capacity(1073741824, x))
        .map(fastx::Reader::new)?;
    let mut writer = stdout();

    for result in reader.records() {
        let mut record = result?;

        record.fastq2fasta();

        record.to_writer(&mut writer)?;
    }

    Ok(())
}
