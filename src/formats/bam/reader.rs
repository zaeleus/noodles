use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

use byteorder::{LittleEndian, ReadBytesExt};
use flate2::read::MultiGzDecoder;

use formats::bam::{Cigar, Data, Flag, Quality, Record, Reference, Sequence};

static MAGIC_NUMBER: &[u8] = b"BAM\x01";

type BamHeader = String;

pub struct Reader<R> {
    reader: MultiGzDecoder<R>,
}

impl<R: Read> Reader<R> {
    pub fn open<P>(path: P) -> io::Result<Reader<File>> where P: AsRef<Path> {
        let file = File::open(path)?;
        let decoder = MultiGzDecoder::new(file);
        Ok(Reader::new(decoder))
    }

    pub fn new(reader: MultiGzDecoder<R>) -> Reader<R> {
        Reader { reader }
    }

    pub fn header(&mut self) -> io::Result<BamHeader> {
        let magic = self.read_magic()?;

        if magic != MAGIC_NUMBER {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid BAM header",
            ));
        }

        let l_text = self.read_l_text()?;
        let text = self.read_text(l_text as usize)?;

        Ok(String::from_utf8(text).unwrap())
    }

    pub fn references(&mut self) -> io::Result<References<R>> {
        let n_ref = self.read_n_ref()?;
        Ok(References::new(self, n_ref as usize))
    }

    pub fn records(&mut self) -> Records<R> {
        Records::new(self)
    }

    /// Reads the BAM magic string
    fn read_magic(&mut self) -> io::Result<Vec<u8>> {
        let mut buf = vec![0; 4];
        self.reader.read_exact(&mut buf)?;
        Ok(buf)
    }

    /// Reads the length of the header
    fn read_l_text(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    /// Reads the header text in SAM
    fn read_text(&mut self, l_text: usize) -> io::Result<Vec<u8>> {
        let mut buf = vec![0; l_text];
        self.reader.read_exact(&mut buf)?;
        Ok(buf)
    }

    /// Reads the number of reference sequences
    fn read_n_ref(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    /// Reads the length of the reference name + 1
    fn read_l_name(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    /// Reads the reference sequence name
    fn read_name(&mut self, l_name: usize) -> io::Result<String> {
        let mut buf = vec![0; l_name];
        self.reader.read_exact(&mut buf).unwrap();
        buf.pop();
        Ok(String::from_utf8(buf).unwrap())
    }

    /// Reads the length of the reference sequence
    fn read_l_ref(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    /// Reads the total length of the alignment block - 4 (this field)
    fn read_block_size(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    /// Reads the index in the list of reference sequences
    fn read_ref_id(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    /// Reads the start position of the alignment (0-based)
    fn read_pos(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    /// Reads the length of the read name
    fn read_l_read_name(&mut self) -> io::Result<u8> {
        self.reader.read_u8()
    }

    /// Reads the mapping quality
    fn read_mapq(&mut self) -> io::Result<u8> {
        self.reader.read_u8()
    }

    /// Reads the BAI bin
    fn read_bin(&mut self) -> io::Result<u16> {
        self.reader.read_u16::<LittleEndian>()
    }

    /// Reads the number of operations in CIGAR
    fn read_n_cigar_op(&mut self) -> io::Result<u16> {
        self.reader.read_u16::<LittleEndian>()
    }

    /// Reads bitwise flags
    fn read_flag(&mut self) -> io::Result<u16> {
        self.reader.read_u16::<LittleEndian>()
    }

    /// Reads the length of the sequence
    fn read_l_seq(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    /// Reads the ref ID of the next sequence
    fn read_next_ref_id(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    /// Reads the start position of the next sequence (0-based)
    fn read_next_pos(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    /// Reads the template length
    fn read_tlen(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    /// Reads the raw read name
    fn read_read_name(&mut self, l_read_name: usize) -> io::Result<Vec<u8>> {
        let mut buf = vec![0; l_read_name];
        self.reader.read_exact(&mut buf)?;
        buf.pop();
        Ok(buf)
    }

    /// Reads the raw CIGAR
    fn read_cigar(&mut self, n_cigar_op: usize) -> io::Result<Vec<u32>> {
        let mut buf = vec![0; n_cigar_op];
        self.reader.read_u32_into::<LittleEndian>(&mut buf)?;
        Ok(buf)
    }

    /// Reads the raw sequence
    fn read_seq(&mut self, l_seq: usize) -> io::Result<Vec<u8>> {
        let len = (l_seq + 1) / 2;
        let mut buf = vec![0; len];
        self.reader.read_exact(&mut buf)?;
        Ok(buf)
    }

    /// Reads the raw quality
    fn read_qual(&mut self, l_seq: usize) -> io::Result<Vec<u8>> {
        let mut buf = vec![0; l_seq];
        self.reader.read_exact(&mut buf)?;
        Ok(buf)
    }

    /// Reads the raw auxiliary data
    fn read_data(&mut self, len: usize) -> io::Result<Vec<u8>> {
        let mut buf = vec![0; len as usize];
        self.reader.read_exact(&mut buf).unwrap();
        Ok(buf)
    }
}

pub struct References<'a, R: 'a> {
    reader: &'a mut Reader<R>,
    i: usize,
    len: usize,
}

impl<'a, R: 'a + Read> References<'a, R> {
    fn new(reader: &'a mut Reader<R>, len: usize) -> References<R> {
        References { reader, i: 0, len }
    }
}

impl<'a, R: 'a + Read> Iterator for References<'a, R> {
    type Item = Reference;

    fn next(&mut self) -> Option<Reference> {
        if self.i >= self.len {
            return None
        }

        let l_name = self.reader.read_l_name().unwrap();
        let name = self.reader.read_name(l_name as usize).unwrap();
        let l_ref = self.reader.read_l_ref().unwrap();

        self.i += 1;

        Some(Reference::new(name, l_ref))
    }
}

pub struct Records<'a, R: 'a> {
    reader: &'a mut Reader<R>,
}

impl<'a, R: 'a + Read> Records<'a, R> {
    fn new(reader: &'a mut Reader<R>) -> Records<R> {
        Records { reader }
    }
}

impl<'a, R: 'a + Read> Iterator for Records<'a, R> {
    type Item = Record;

    fn next(&mut self) -> Option<Record> {
        let block_size = match self.reader.read_block_size() {
            Ok(bs) => bs,
            Err(_) => return None,
        };

        let ref_id = self.reader.read_ref_id().unwrap();
        let pos = self.reader.read_pos().unwrap();
        let l_read_name = self.reader.read_l_read_name().unwrap();
        let mapq = self.reader.read_mapq().unwrap();
        let _bin = self.reader.read_bin().unwrap();
        let n_cigar_op = self.reader.read_n_cigar_op().unwrap();
        let flag = self.reader.read_flag().unwrap();
        let l_seq = self.reader.read_l_seq().unwrap();
        let next_ref_id = self.reader.read_next_ref_id().unwrap();
        let next_pos = self.reader.read_next_pos().unwrap();
        let _tlen = self.reader.read_tlen().unwrap();

        let read_name = self.reader
            .read_read_name(l_read_name as usize)
            .map(|n| String::from_utf8(n).unwrap())
            .unwrap();

        let cigar = self.reader.read_cigar(n_cigar_op as usize).unwrap();
        let seq = self.reader.read_seq(l_seq as usize).unwrap();
        let qual = self.reader.read_qual(l_seq as usize).unwrap();

        // aux data
        let bytes_read = 32
            + l_read_name as i32
            + 4 * n_cigar_op as i32
            + (l_seq as i32 + 1) / 2
            + l_seq;

        let len = block_size - bytes_read;
        let data = self.reader.read_data(len as usize).unwrap();

        Some(Record::new(
            ref_id,
            pos,
            mapq,
            Flag::new(flag),
            next_ref_id,
            next_pos,
            read_name,
            Cigar::new(cigar),
            Sequence::new(seq, l_seq as usize),
            Quality::new(qual),
            Data::new(data),
        ))
    }
}
