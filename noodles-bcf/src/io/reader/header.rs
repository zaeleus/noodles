mod format_version;
pub(crate) mod magic_number;

use std::io::{self, BufRead, BufReader, Read, Take};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::{self as vcf, header::StringMaps};

pub(super) use self::{format_version::read_format_version, magic_number::read_magic_number};

#[derive(Clone, Copy, Debug)]
enum State {
    Length,
    Header,
    Padding,
    Done,
}

struct Reader<R> {
    inner: BufReader<Take<R>>,
    state: State,
}

impl<R> Reader<R>
where
    R: Read,
{
    fn new(inner: R) -> Self {
        Self {
            inner: BufReader::new(inner.take(4)),
            state: State::Length,
        }
    }
}

impl<R> Read for Reader<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut src = self.fill_buf()?;
        let amt = src.read(buf)?;
        self.consume(amt);
        Ok(amt)
    }
}

impl<R> BufRead for Reader<R>
where
    R: Read,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        const NUL: u8 = 0x00;

        loop {
            match self.state {
                State::Length => {
                    let header_len = self.inner.read_u32::<LittleEndian>().map(u64::from)?;
                    self.inner.get_mut().set_limit(header_len);
                    self.state = State::Header;
                }
                State::Header => {
                    let src = self.inner.fill_buf()?;

                    if src.is_empty() {
                        self.state = State::Done;
                    } else if src[0] == NUL {
                        self.state = State::Padding;
                    } else if let Some(i) = src.iter().position(|&b| b == NUL) {
                        // NOTE: This is a borrowck workaround. See NLL problem case #3.
                        let src = self.inner.fill_buf()?;
                        return Ok(&src[..i]);
                    } else {
                        // NOTE: This is a borrowck workaround. See NLL problem case #3.
                        let src = self.inner.fill_buf()?;
                        return Ok(src);
                    }
                }
                State::Padding => {
                    io::copy(&mut self.inner, &mut io::sink())?;
                    self.state = State::Done;
                }
                State::Done => return Ok(&[]),
            }
        }
    }

    fn consume(&mut self, amt: usize) {
        self.inner.consume(amt);
    }
}

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<vcf::Header>
where
    R: Read,
{
    let mut parser = vcf::header::Parser::default();
    let mut string_maps = StringMaps::default();

    let mut reader = Reader::new(reader);
    let mut buf = Vec::new();

    while read_line(&mut reader, &mut buf)? != 0 {
        let entry = parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        string_maps
            .insert_entry(&entry)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    let mut header = parser
        .finish()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *header.string_maps_mut() = string_maps;

    Ok(header)
}

fn read_line<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    dst.clear();

    match reader.read_until(LINE_FEED, dst)? {
        0 => Ok(0),
        n => {
            if dst.ends_with(&[LINE_FEED]) {
                dst.pop();

                if dst.ends_with(&[CARRIAGE_RETURN]) {
                    dst.pop();
                }
            }

            Ok(n)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_header() -> io::Result<()> {
        use vcf::header::FileFormat;

        const NUL: u8 = 0x00;

        let raw_header = b"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
";

        let mut data = 61u32.to_le_bytes().to_vec(); // l_text
        data.extend_from_slice(raw_header);
        data.push(NUL);

        let mut reader = &data[..];
        let actual = read_header(&mut reader)?;

        let expected = vcf::Header::builder()
            .set_file_format(FileFormat::new(4, 3))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
