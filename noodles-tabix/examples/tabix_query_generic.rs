//! Queries an indexed TSV with the given region.
//!
//! The input must have an associated tabix index in the same directory.
//!
//! The result matches the output of `tabix <src> <region>`.

use std::{env, fs::File};

use noodles_csi as csi;
use noodles_tabix as tabix;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");
    let region = args.next().expect("missing region").parse()?;

    let tabix_src = format!("{src}.tbi");
    let index = tabix::read(tabix_src)?;
    let mut reader = File::open(src).map(|f| csi::io::IndexedReader::new(f, index))?;

    let query = reader.query(&region)?;

    for result in query {
        let record = result?;
        println!("{}", record.as_ref());
    }

    Ok(())
}
