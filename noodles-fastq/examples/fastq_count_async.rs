//! Counts the number of records in a FASTQ file.

use std::env;

use noodles_fastq as fastq;
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
        .map(fastq::AsyncReader::new)?;

    let mut record = fastq::Record::default();
    let mut n = 0;

    loop {
        match reader.read_record(&mut record).await? {
            0 => break,
            _ => n += 1,
        }
    }

    println!("{}", n);

    Ok(())
}
