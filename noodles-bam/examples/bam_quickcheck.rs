//! Performs very simple validations on a BAM file.
//!
//! The checks are similar to the behavior of `samtools quickcheck <src>`.

use std::{
    env,
    fs::File,
    io::{self, Read, Seek, SeekFrom},
};

use noodles_bam as bam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;

    if header.reference_sequences().is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "empty reference sequence dictionary",
        ));
    }

    let mut file = reader.into_inner().into_inner();

    if !has_eof(&mut file)? {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing EOF marker",
        ));
    }

    Ok(())
}

fn has_eof<R>(reader: &mut R) -> io::Result<bool>
where
    R: Read + Seek,
{
    // _Sequence Alignment/Map Format Specification_ (2024-11-06) § 4.1.2 "End-of-file marker".
    const EOF: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    const EOF_OFFSET: SeekFrom = SeekFrom::End(-(EOF.len() as i64));

    reader.seek(EOF_OFFSET)?;

    let mut buf = [0; EOF.len()];
    reader.read_exact(&mut buf)?;

    Ok(buf == EOF)
}
