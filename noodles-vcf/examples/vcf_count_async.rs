//! Counts the number of records in a VCF file.
//!
//! The result matches the output of `bcftools view --no-header <src> | wc -l`.

use std::env;

use noodles_vcf as vcf;
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
        .map(vcf::AsyncReader::new)?;

    reader.read_header().await?;

    let mut buf = String::new();
    let mut n = 0;

    loop {
        match reader.read_record(&mut buf).await {
            Ok(0) => break,
            Ok(_) => n += 1,
            Err(e) => return Err(e),
        }
    }

    println!("{}", n);

    Ok(())
}
