use std::{
    cmp,
    io::{self, Write},
};

use byteorder::{LittleEndian, WriteBytesExt};

use flate2::{write::DeflateEncoder, Compression, Crc};

const MAX_BGZF_BLOCK_SIZE: u32 = 65536; // bytes

#[non_exhaustive]
enum CompressionMethod {
    Deflate = 8,
}

#[non_exhaustive]
enum OperatingSystem {
    Unknown = 255,
}

const GZIP_ID1: u8 = 0x1f;
const GZIP_ID2: u8 = 0x8b;
const GZIP_FLG_FEXTRA: u8 = 0x04;
const GZIP_MTIME_NONE: u32 = 0;
const GZIP_XFL_NONE: u8 = 0x00;
const GZIP_XLEN_BGZF: u16 = 6;

const BGZF_SI1: u8 = 0x42;
const BGZF_SI2: u8 = 0x43;
const BGZF_SLEN: u16 = 2;

#[derive(Debug)]
pub struct Writer<W>
where
    W: Write,
{
    inner: W,
    encoder: DeflateEncoder<Vec<u8>>,
    crc: Crc,
}

impl<W> Writer<W>
where
    W: Write,
{
    pub fn new(inner: W) -> Self {
        Self {
            inner,
            encoder: DeflateEncoder::new(Vec::new(), Compression::default()),
            crc: Crc::new(),
        }
    }

    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    fn flush_block(&mut self) -> io::Result<()> {
        self.encoder.try_finish()?;
        let data = self.encoder.get_ref();

        write_header(&mut self.inner, data.len() as u64)?;
        self.inner.write_all(&data[..])?;
        write_trailer(&mut self.inner, self.crc.sum(), self.crc.amount())?;

        self.encoder.reset(Vec::new())?;
        self.crc.reset();

        Ok(())
    }
}

impl<W> Write for Writer<W>
where
    W: Write,
{
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let total_uncompressed_bytes_written = self.crc.amount();

        if total_uncompressed_bytes_written >= MAX_BGZF_BLOCK_SIZE {
            self.flush()?;
            return Err(io::Error::from(io::ErrorKind::Interrupted));
        }

        let bytes_to_be_written = cmp::min(
            (MAX_BGZF_BLOCK_SIZE - total_uncompressed_bytes_written) as usize,
            buf.len(),
        );
        let bytes_written = self.encoder.write(&buf[..bytes_to_be_written])?;
        self.crc.update(&buf[..bytes_written]);

        Ok(bytes_written)
    }

    fn flush(&mut self) -> io::Result<()> {
        if self.crc.amount() > 0 {
            self.flush_block()
        } else {
            Ok(())
        }
    }
}

impl<W> Drop for Writer<W>
where
    W: Write,
{
    fn drop(&mut self) {
        // Ignore a failed flush.
        //
        // Interestingly, this matches the behavior of `std::io::BufWriter`.
        let _r = self.flush();
    }
}

pub fn write_header<W>(writer: &mut W, block_size: u64) -> io::Result<()>
where
    W: Write,
{
    writer.write_u8(GZIP_ID1)?;
    writer.write_u8(GZIP_ID2)?;
    writer.write_u8(CompressionMethod::Deflate as u8)?;
    writer.write_u8(GZIP_FLG_FEXTRA)?;
    writer.write_u32::<LittleEndian>(GZIP_MTIME_NONE)?;
    writer.write_u8(GZIP_XFL_NONE)?;
    writer.write_u8(OperatingSystem::Unknown as u8)?;
    writer.write_u16::<LittleEndian>(GZIP_XLEN_BGZF)?;

    writer.write_u8(BGZF_SI1)?;
    writer.write_u8(BGZF_SI2)?;
    writer.write_u16::<LittleEndian>(BGZF_SLEN)?;

    let bsize = (block_size as u16) + 18 + 8 - 1;
    writer.write_u16::<LittleEndian>(bsize)?;

    Ok(())
}

pub fn write_trailer<W>(writer: &mut W, checksum: u32, uncompressed_size: u32) -> io::Result<()>
where
    W: Write,
{
    writer.write_u32::<LittleEndian>(checksum)?;
    writer.write_u32::<LittleEndian>(uncompressed_size)?;
    Ok(())
}
