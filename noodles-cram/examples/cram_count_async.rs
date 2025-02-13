//! Counts the number of records in a CRAM file.
//!
//! Note that this example counts each successfully read record and not by summing the number of
//! records field in each container or slice header.
//!
//! The result matches the output of `samtools view --count <src>`.

use std::env;

use noodles_cram::{self as cram, io::reader::Container};
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(cram::r#async::io::Reader::new)?;
    reader.read_header().await?;

    let mut container = Container::default();
    let mut n = 0;

    while reader.read_container(&mut container).await? != 0 {
        let compression_header = container.compression_header()?;

        for result in container.slices() {
            let slice = result?;

            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

            let records =
                slice.records(&compression_header, &core_data_src, &external_data_srcs)?;

            n += records.len();
        }
    }

    println!("{n}");

    Ok(())
}
