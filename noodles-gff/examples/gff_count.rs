//! Counts the number of records in a GFF file.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_gff as gff;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(BufReader::new).map(gff::Reader::new)?;
    reader.read_header()?;

    let mut n = 0;
    let mut buf = String::new();

    loop {
        buf.clear();

        match reader.read_record(&mut buf) {
            Ok(0) => break,
            Ok(_) => n += 1,
            Err(e) => return Err(e),
        }
    }

    println!("{}", n);

    Ok(())
}
