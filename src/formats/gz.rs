use std::io::{self, BufRead, Read};

use byteorder::{LittleEndian, ReadBytesExt};
use flate2::bufread::DeflateDecoder;

static GZIP_MAGIC_NUMBER: &[u8] = &[0x1f, 0x8b];

const DEFAULT_COMPRESSION_METHOD: u8 = 8; // DEFLATE
const FHCRC: u8 = 1 << 1;
const FEXTRA: u8 = 1 << 2;
const FNAME: u8 = 1 << 3;
const FCOMMENT: u8 = 1 << 4;

pub struct MultiGzDecoder<R> {
    inner: DeflateDecoder<R>,
    header: io::Result<()>,
    finished: bool,
}

impl<R: BufRead> MultiGzDecoder<R> {
    pub fn new(mut reader: R) -> MultiGzDecoder<R> {
        let header = read_gz_header(&mut reader);
        let inner = DeflateDecoder::new(reader);

        MultiGzDecoder {
            inner,
            header,
            finished: false,
        }
    }

    fn finish_member(&mut self) -> io::Result<usize> {
        if self.finished {
            return Ok(0);
        }

        let mut trailer = [0; 8];
        self.inner.get_mut().read_exact(&mut trailer)?;

        let remaining = {
            let buf = self.inner.get_mut().fill_buf()?;

            if buf.is_empty() {
                self.finished = true;
                return Ok(0);
            } else {
                buf.len()
            }
        };

        self.header = read_gz_header(self.inner.get_mut());

        self.inner.reset_data();

        Ok(remaining)
    }
}

impl<R: BufRead> Read for MultiGzDecoder<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self.inner.read(buf)? {
            0 => match self.finish_member() {
                Ok(0) => Ok(0),
                Ok(_) => self.read(buf),
                Err(e) => Err(e),
            },
            n => Ok(n),
        }
    }
}

fn read_gz_header<R: BufRead>(reader: &mut R) -> io::Result<()> {
    let mut id = [0; 2];
    reader.read_exact(&mut id)?;
    if id != GZIP_MAGIC_NUMBER {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            String::from("gz: invalid gzip header")
        ));
    }

    let cm = reader.read_u8()?;
    if cm != DEFAULT_COMPRESSION_METHOD {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            String::from("gz: unsupported compression method"),
        ));
    }

    let flg = reader.read_u8()?;

    let _mtime = reader.read_u32::<LittleEndian>()?;
    let _xfl = reader.read_u8()?;
    let _os = reader.read_u8()?;

    if flg & FEXTRA != 0 {
        let xlen = reader.read_u16::<LittleEndian>()?;
        let mut buf = vec![0; xlen as usize];
        reader.read_exact(&mut buf)?;
    }

    if flg & FNAME != 0 {
        let mut buf = Vec::new();
        reader.read_until(b'\0', &mut buf)?;
    }

    if flg & FCOMMENT != 0 {
        let mut buf = Vec::new();
        reader.read_until(b'\0', &mut buf)?;
    }

    if flg & FHCRC != 0 {
        reader.read_u16::<LittleEndian>()?;
    }

    Ok(())
}
