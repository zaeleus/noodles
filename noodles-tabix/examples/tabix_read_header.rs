//! Prints the header of the file associated with the index.
//!
//! The results match the output of `tabix --only-header <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufRead},
};

use noodles_bgzf as bgzf;
use noodles_csi::BinningIndex;
use noodles_tabix as tabix;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let tabix_src = format!("{src}.tbi");
    let index = tabix::fs::read(tabix_src)?;

    let header = index
        .header()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing tabix header"))?;

    let reader = File::open(src).map(bgzf::io::Reader::new)?;
    let line_comment_prefix = char::from(header.line_comment_prefix());

    for result in reader.lines() {
        let line = result?;

        if !line.starts_with(line_comment_prefix) {
            break;
        }

        println!("{line}");
    }

    Ok(())
}
