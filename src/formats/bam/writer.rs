use std::{
    fs::File,
    io::{self, Write},
    path::Path,
};

use byteorder::{LittleEndian, WriteBytesExt};

use crate::formats::bgzf::BgzfEncoder;

use super::{Record, Reference, MAGIC_NUMBER};

pub struct Writer<W: Write> {
    writer: BgzfEncoder<W>,
}

impl<W: Write> Writer<W> {
    pub fn create<P>(dst: P) -> io::Result<Writer<File>>
    where
        P: AsRef<Path>,
    {
        let encoder = BgzfEncoder::<File>::create(dst)?;
        Ok(Writer::new(encoder))
    }

    pub fn new(writer: BgzfEncoder<W>) -> Writer<W> {
        Writer { writer }
    }

    pub fn write_header(&mut self, text: &str) -> io::Result<()> {
        self.writer.write_all(MAGIC_NUMBER)?;

        self.writer.write_i32::<LittleEndian>(text.len() as i32)?;
        self.writer.write_all(text.as_bytes())?;

        Ok(())
    }

    pub fn write_footer(&mut self) -> io::Result<()> {
        self.writer.write_footer()
    }

    pub fn write_references(&mut self, references: &[Reference]) -> io::Result<()> {
        self.writer
            .write_i32::<LittleEndian>(references.len() as i32)?;

        for reference in references {
            let len = reference.name().len() + 1;
            self.writer.write_i32::<LittleEndian>(len as i32)?;

            self.writer.write_all(reference.name().as_bytes())?;
            self.writer.write_all(b"\0")?;

            self.writer.write_i32::<LittleEndian>(reference.len())?;
        }

        Ok(())
    }

    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        self.writer.write_i32::<LittleEndian>(record.block_size())?;

        self.writer.write_i32::<LittleEndian>(record.ref_id())?;
        self.writer.write_i32::<LittleEndian>(record.pos())?;
        self.writer.write_u8(record.l_read_name())?;
        self.writer.write_u8(record.mapq())?;
        self.writer.write_u16::<LittleEndian>(record.bin())?;
        self.writer.write_u16::<LittleEndian>(record.n_cigar_op())?;
        self.writer
            .write_u16::<LittleEndian>(record.flag().inner())?;
        self.writer.write_i32::<LittleEndian>(record.l_seq())?;
        self.writer
            .write_i32::<LittleEndian>(record.next_ref_id())?;
        self.writer.write_i32::<LittleEndian>(record.next_pos())?;
        self.writer.write_i32::<LittleEndian>(record.tlen())?;

        self.writer.write_all(record.read_name())?;
        self.writer.write_all(b"\0")?;

        self.writer.write_all(&record.cigar())?;
        self.writer.write_all(&record.seq())?;
        self.writer.write_all(&record.qual())?;
        self.writer.write_all(&record.data())?;

        Ok(())
    }

    pub fn finish_block(&mut self) -> io::Result<()> {
        self.writer.finish_block()
    }
}

impl<W: Write> Write for Writer<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.writer.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}
