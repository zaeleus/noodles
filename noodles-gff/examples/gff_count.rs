//! Counts the number of records in a GFF file.

use std::{env, fs::File, io::BufReader};

use noodles_gff as gff;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(BufReader::new).map(gff::Reader::new)?;

    let mut n = 0;
    let mut buf = String::new();

    loop {
        buf.clear();

        match reader.read_line(&mut buf) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => return Err(e.into()),
        }

        let line: gff::Line = buf.parse()?;

        match line {
            gff::Line::Record(_) => n += 1,
            gff::Line::Directive(d) => {
                if let gff::Directive::StartOfFasta = d {
                    break;
                }
            }
            _ => {}
        }
    }

    println!("{}", n);

    Ok(())
}
