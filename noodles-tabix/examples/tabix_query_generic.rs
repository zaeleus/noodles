//! Queries an indexed TSV with the given region.
//!
//! The input must have an associated tabix index in the same directory.
//!
//! The result matches the output of `tabix <src> <region>`.

use std::env;

use noodles_tabix as tabix;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");
    let region = args.next().expect("missing region").parse()?;

    let mut reader = tabix::io::indexed_reader::Builder::default().build_from_path(src)?;
    let query = reader.query(&region)?;

    for result in query {
        let record = result?;
        println!("{}", record.as_ref());
    }

    Ok(())
}
