//! Variant format filesystem operations.

mod file;

use std::{
    io::{self, BufRead},
    path::Path,
};

use noodles_bcf as bcf;
use noodles_vcf as vcf;

use self::file::File;
use super::io::{Format, Reader};

/// Opens a variant format file.
///
/// # Examples
///
/// ```no_run
/// use noodles_util::variant;
/// let reader = variant::fs::open("in.vcf")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn open<P>(src: P) -> io::Result<Reader<File>>
where
    P: AsRef<Path>,
{
    let mut file = File::open(src)?;

    let inner = match detect_format(&mut file)? {
        Format::Vcf => super::io::reader::Inner::Vcf(vcf::io::Reader::new(file)),
        Format::Bcf => super::io::reader::Inner::Bcf(bcf::io::Reader::from(file)),
    };

    Ok(Reader { inner })
}

fn detect_format<R>(reader: &mut R) -> io::Result<Format>
where
    R: BufRead,
{
    const BCF_MAGIC_NUMBER: [u8; 3] = *b"BCF";

    let src = reader.fill_buf()?;

    if let Some(buf) = src.get(..BCF_MAGIC_NUMBER.len()) {
        if buf == BCF_MAGIC_NUMBER {
            return Ok(Format::Bcf);
        }
    }

    Ok(Format::Vcf)
}
