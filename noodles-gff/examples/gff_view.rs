//! Prints all lines in a GFF file, up to the FASTA section.
//!
//! Lines are parsed as either a directive, comment, or record.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_gff::{self as gff, Directive, Line};

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(BufReader::new).map(gff::Reader::new)?;

    let stdout = io::stdout().lock();
    let mut writer = gff::Writer::new(stdout);

    for result in reader.lines() {
        let line = result?;

        if line == Line::Directive(Directive::StartOfFasta) {
            break;
        }

        writer.write_line(&line)?;
    }

    Ok(())
}
