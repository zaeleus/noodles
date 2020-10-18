use structopt::StructOpt;
use std::{fs::File, io::BufReader};
/// A tutorial on noodles bioinformatics I/O
#[derive(StructOpt)]
struct Cli {
    /// Path to the file to read
    #[structopt(
    short="i", 
    long = "--infile",
    parse(from_os_str)
)]
    infile: std::path::PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::from_args();
    let mut reader = File::open(&args.infile).map(BufReader::new).map(noodles_sam::Reader::new)?;
    for result in reader.records() {
        let record = result;
        println!("{:?}", record);
 }
    Ok(())
}