use std::io::{self, Write};
use std::fs::File;
use std::path::Path;

use byteorder::{LittleEndian, WriteBytesExt};
use flate2::{Crc, Compression};
use flate2::write::DeflateEncoder;

const DEFAULT_COMPRESSION_METHOD: u8 = 8; // DEFLATE
const FEXTRA: u8 = 1 << 2;
const DEFAULT_FLAGS: u8 = FEXTRA;
const DEFAULT_MODIFICATION_TIME: u32 = 0;
const DEFAULT_EXTRA_FLAGS: u8 = 0;
const DEFAULT_OPERATING_SYSTEM: u8 = 255; // unknown
const BGZF_SUBFIELD_LENGTH: u16 = 2; // bytes

const MAX_BLOCK_SIZE: u64 = 65536; // bytes

static GZIP_MAGIC_NUMBER: &[u8] = &[0x1f, 0x8b];
static BGZF_SUBFIELD_ID: &[u8] = b"BC";

static EOF: &[u8] = &[
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
];

#[derive(Debug)]
pub struct BgzfEncoder<W: Write> {
    crc: Crc,
    encoder: DeflateEncoder<Vec<u8>>,
    writer: W,
}

impl<W: Write> BgzfEncoder<W> {
    pub fn create<P>(dst: P) -> io::Result<BgzfEncoder<File>> where P: AsRef<Path> {
        let file = File::create(dst)?;
        Ok(BgzfEncoder::new(file))
    }

    pub fn new(writer: W) -> BgzfEncoder<W> {
        let crc = Crc::new();

        let buf = Vec::with_capacity(MAX_BLOCK_SIZE as usize);
        let encoder = DeflateEncoder::new(buf, Compression::default());

        BgzfEncoder { crc, encoder, writer }
    }

    fn block_header(&self, block_size: u16) -> Vec<u8> {
        assert!(block_size > 0);

        let mut header = Vec::new();

        header.write_all(GZIP_MAGIC_NUMBER).unwrap();
        header.write_u8(DEFAULT_COMPRESSION_METHOD).unwrap();
        header.write_u8(DEFAULT_FLAGS).unwrap();
        header.write_u32::<LittleEndian>(DEFAULT_MODIFICATION_TIME).unwrap();
        header.write_u8(DEFAULT_EXTRA_FLAGS).unwrap();
        header.write_u8(DEFAULT_OPERATING_SYSTEM).unwrap();

        let mut extra = Vec::new();

        extra.write_all(BGZF_SUBFIELD_ID).unwrap();
        extra.write_u16::<LittleEndian>(BGZF_SUBFIELD_LENGTH).unwrap();
        extra.write_u16::<LittleEndian>(block_size - 1).unwrap();

        let xlen = extra.len() as u16;
        header.write_u16::<LittleEndian>(xlen).unwrap();

        header.extend(extra);

        header
    }

    pub fn finish_block(&mut self) -> io::Result<()> {
        if self.encoder.total_in() == 0 {
            return Ok(());
        }

        self.encoder.try_finish().unwrap();

        {
            let data = self.encoder.get_ref();

            let len = data.len() as u16;
            let header = self.block_header(len + 6 + 20);
            self.writer.write_all(&header)?;

            self.writer.write_all(data)?;
        }

        let checksum = self.crc.sum();
        self.writer.write_u32::<LittleEndian>(checksum)?;

        let input_size = self.encoder.total_in() as u32;
        self.writer.write_u32::<LittleEndian>(input_size)?;

        self.reset()
    }

    pub fn write_footer(&mut self) -> io::Result<()> {
        self.finish_block()?;
        self.writer.write_all(EOF)
    }

    fn reset(&mut self) -> io::Result<()> {
        self.crc = Crc::new();
        let buf = Vec::with_capacity(MAX_BLOCK_SIZE as usize);
        self.encoder.reset(buf)?;
        Ok(())
    }
}

impl<W: Write> Write for BgzfEncoder<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let len = buf.len() as u64;
        let next_size = len + self.encoder.total_in();

        if next_size >= MAX_BLOCK_SIZE {
            self.finish_block()?;
        }

        let bytes_written = self.encoder.write(buf)?;
        self.crc.update(&buf[..bytes_written]);
        Ok(bytes_written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.finish_block()?;
        self.writer.flush()
    }
}

impl<W: Write> Drop for BgzfEncoder<W> {
    fn drop(&mut self) {
        self.flush().unwrap();
    }
}
