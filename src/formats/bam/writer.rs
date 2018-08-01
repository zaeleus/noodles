use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use byteorder::{LittleEndian, WriteBytesExt};

use formats::bam::MAGIC_NUMBER;
use formats::bam::{Record, Reference};

pub struct Writer<W> {
    writer: W,
}

impl<W: Write> Writer<W> {
    pub fn create<P>(dst: P) -> io::Result<Writer<BufWriter<File>>>
    where
        P: AsRef<Path>,
    {
        let file = File::create(dst)?;
        let writer = BufWriter::new(file);
        Ok(Writer::new(writer))
    }

    pub fn new(writer: W) -> Writer<W> {
        Writer { writer }
    }

    pub fn write_header(&mut self, text: &str) -> io::Result<()> {
        self.writer.write_all(MAGIC_NUMBER)?;

        self.writer.write_i32::<LittleEndian>(text.len() as i32)?;
        self.writer.write_all(text.as_bytes())?;

        Ok(())
    }

    pub fn write_references(&mut self, references: &[Reference]) -> io::Result<()> {
        self.writer.write_i32::<LittleEndian>(references.len() as i32)?;

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
        self.writer.write_u16::<LittleEndian>(record.flag().inner())?;
        self.writer.write_i32::<LittleEndian>(record.l_seq())?;
        self.writer.write_i32::<LittleEndian>(record.next_ref_id())?;
        self.writer.write_i32::<LittleEndian>(record.next_ref_pos())?;
        self.writer.write_i32::<LittleEndian>(record.tlen())?;

        self.writer.write_all(record.read_name().as_bytes())?;
        self.writer.write_all(b"\0")?;

        for &u in record.cigar().iter() {
            self.writer.write_u32::<LittleEndian>(u)?;
        }

        self.writer.write_all(record.sequence().as_bytes())?;
        self.writer.write_all(record.quality().as_bytes())?;
        self.writer.write_all(record.data().as_bytes())?;

        Ok(())
    }
}
