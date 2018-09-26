use std::fs::File;
use std::io::{self, BufReader, Read};
use std::path::Path;

use byteorder::{LittleEndian, ReadBytesExt};

use formats::bam::MAGIC_NUMBER;
use formats::bam::{ByteRecord, Record, Reference};
use formats::gz::MultiGzDecoder;

type BamHeader = Vec<u8>;

pub struct Reader<R> {
    reader: MultiGzDecoder<BufReader<R>>,
}

impl<R: Read> Reader<R> {
    pub fn open<P>(path: P) -> io::Result<Reader<File>> where P: AsRef<Path> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let decoder = MultiGzDecoder::new(reader);
        Ok(Reader::new(decoder))
    }

    pub fn new(reader: MultiGzDecoder<BufReader<R>>) -> Reader<R> {
        Reader { reader }
    }

    pub fn header(&mut self) -> io::Result<BamHeader> {
        let mut magic = [0; 4];
        self.reader.read_exact(&mut magic)?;

        if magic != MAGIC_NUMBER {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid BAM header",
            ));
        }

        let l_text = self.read_l_text()?;

        let mut text = vec![0; l_text as usize];
        self.reader.read_exact(&mut text)?;

        Ok(text)
    }

    pub fn references(&mut self) -> io::Result<References<R>> {
        let n_ref = self.read_n_ref()?;
        Ok(References::new(self, n_ref as usize))
    }

    pub fn records(&mut self) -> Records<R> {
        Records::new(self)
    }

    fn read_reference(&mut self) -> io::Result<Reference> {
        let l_name = self.read_l_name()?;
        let name = self.read_name(l_name as usize)?;
        let l_ref = self.read_l_ref()?;
        Ok(Reference::new(name, l_ref))
    }

    pub fn read_byte_record(&mut self, record: &mut ByteRecord) -> io::Result<usize> {
        let block_size = match self.read_block_size() {
            Ok(bs) => bs as usize,
            Err(_) => return Ok(0),
        };

        record.resize(block_size);
        self.reader.read_exact(record)?;

        Ok(block_size)
    }

    fn read_l_text(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    fn read_n_ref(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    fn read_l_name(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    fn read_name(&mut self, l_name: usize) -> io::Result<String> {
        let mut buf = vec![0; l_name];
        self.reader.read_exact(&mut buf)?;
        buf.pop();
        Ok(String::from_utf8(buf).unwrap())
    }

    fn read_l_ref(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    fn read_block_size(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
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
    type Item = io::Result<Reference>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i >= self.len {
            return None;
        }

        let result = self.reader.read_reference();
        self.i += 1;
        Some(result)
    }
}

pub struct Records<'a, R: 'a> {
    reader: &'a mut Reader<R>,
    buf: ByteRecord,
}

impl<'a, R: 'a + Read> Records<'a, R> {
    fn new(reader: &'a mut Reader<R>) -> Records<R> {
        Records {
            reader,
            buf: ByteRecord::new(),
        }
    }
}

impl<'a, R: 'a + Read> Iterator for Records<'a, R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_byte_record(&mut self.buf) {
            Ok(0) => None,
            Ok(_) => Some(Ok(Record::from(&self.buf))),
            Err(e) => Some(Err(e)),
        }
    }
}
