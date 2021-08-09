//! Prints the header of a BCF file.

use std::env;

use noodles_bcf as bcf;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(bcf::AsyncReader::new)?;
    reader.read_file_format().await?;

    let header = reader.read_header().await?;
    print!("{}", header);

    Ok(())
}
