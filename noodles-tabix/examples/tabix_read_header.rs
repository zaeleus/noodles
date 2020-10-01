//! Prints the header of the file associated with the index.
//!
//! The results match the output of `tabix --only-header <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufRead, BufReader},
};

use noodles_bgzf as bgzf;
use noodles_tabix as tabix;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let tabix_src = format!("{}.tbi", src);
    let index = tabix::read(tabix_src)?;

    let reader = File::open(src).map(bgzf::Reader::new).map(BufReader::new)?;
    let line_comment_prefix = char::from(index.header().line_comment_prefix());

    for result in reader.lines() {
        let line = result?;

        if line.starts_with(line_comment_prefix) {
            println!("{}", line);
        } else {
            break;
        }
    }

    Ok(())
}
