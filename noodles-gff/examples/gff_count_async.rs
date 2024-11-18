//! Counts the number of records in a GFF file.
//!
//! Assuming the input does not have a FASTA section, the result matches the output of `grep
//! --count --invert-match '^#' <src>`.

use std::env;

use noodles_gff as gff;
use tokio::{
    fs::File,
    io::{self, BufReader},
};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .await
        .map(BufReader::new)
        .map(gff::r#async::io::Reader::new)?;

    let mut line = gff::Line::default();
    let mut n = 0;

    while reader.read_line(&mut line).await? != 0 {
        if line.as_record().is_some() {
            n += 1;
        }
    }

    println!("{n}");

    Ok(())
}
