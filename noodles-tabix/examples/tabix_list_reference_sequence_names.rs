//! Prints the reference sequence names stored in the index.
//!
//! The results match the output of `tabix --list-chroms <src>`.

use std::{env, io};

use noodles_tabix as tabix;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let tabix_src = format!("{}.tbi", src);
    let index = tabix::read(tabix_src)?;

    for reference_sequence_name in index.reference_sequence_names() {
        println!("{}", reference_sequence_name);
    }

    Ok(())
}
