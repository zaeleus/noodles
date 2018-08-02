use std::io::{self, BufRead, Cursor};

use byteorder::{LittleEndian, ReadBytesExt};

struct Reader<R> {
    reader: R,
}

impl<R: BufRead> Reader<R> {
    fn new(reader: R) -> Reader<R> {
        Reader { reader }
    }

    fn fields(&mut self) -> Fields<R> {
        Fields { reader: self }
    }

    fn read_tag(&mut self) -> io::Result<Vec<u8>> {
        let mut buf = vec![0; 2];
        self.reader.read_exact(&mut buf)?;
        Ok(buf)
    }

    fn read_char(&mut self) -> io::Result<char> {
        self.reader.read_u8().map(|b| b as char)
    }

    fn read_i8(&mut self) -> io::Result<i8> {
        self.reader.read_i8()
    }

    fn read_u8(&mut self) -> io::Result<u8> {
        self.reader.read_u8()
    }

    fn read_i16(&mut self) -> io::Result<i16> {
        self.reader.read_i16::<LittleEndian>()
    }

    fn read_u16(&mut self) -> io::Result<u16> {
        self.reader.read_u16::<LittleEndian>()
    }

    fn read_i32(&mut self) -> io::Result<i32> {
        self.reader.read_i32::<LittleEndian>()
    }

    fn read_u32(&mut self) -> io::Result<u32> {
        self.reader.read_u32::<LittleEndian>()
    }

    fn read_f32(&mut self) -> io::Result<f32> {
        self.reader.read_f32::<LittleEndian>()
    }

    fn read_string(&mut self) -> io::Result<Vec<u8>> {
        let mut buf = Vec::new();
        self.reader.read_until(b'\0', &mut buf)?;
        buf.pop();
        Ok(buf)
    }

    fn read_hex(&mut self) -> io::Result<Vec<u8>> {
        self.read_string()
    }

    fn read_i8_array(&mut self, len: usize) -> io::Result<Vec<i8>> {
        let mut buf = vec![0; len];
        self.reader.read_exact(&mut buf)?;
        Ok(buf.iter().map(|&b| b as i8).collect())
    }

    fn read_u8_array(&mut self, len: usize) -> io::Result<Vec<u8>> {
        let mut buf = vec![0; len];
        self.reader.read_exact(&mut buf)?;
        Ok(buf)
    }

    fn read_i16_array(&mut self, len: usize) -> io::Result<Vec<i16>> {
        let mut buf = vec![0; len];
        self.reader.read_i16_into::<LittleEndian>(&mut buf)?;
        Ok(buf)
    }

    fn read_u16_array(&mut self, len: usize) -> io::Result<Vec<u16>> {
        let mut buf = vec![0; len];
        self.reader.read_u16_into::<LittleEndian>(&mut buf)?;
        Ok(buf)
    }

    fn read_i32_array(&mut self, len: usize) -> io::Result<Vec<i32>> {
        let mut buf = vec![0; len];
        self.reader.read_i32_into::<LittleEndian>(&mut buf)?;
        Ok(buf)
    }

    fn read_u32_array(&mut self, len: usize) -> io::Result<Vec<u32>> {
        let mut buf = vec![0; len];
        self.reader.read_u32_into::<LittleEndian>(&mut buf)?;
        Ok(buf)
    }

    fn read_f32_array(&mut self, len: usize) -> io::Result<Vec<f32>> {
        let mut buf = vec![0.0; len];
        self.reader.read_f32_into::<LittleEndian>(&mut buf)?;
        Ok(buf)
    }
}

struct Fields<'a, R: 'a> {
    reader: &'a mut Reader<R>,
}

impl<'a, R: 'a + BufRead> Iterator for Fields<'a, R> {
    type Item = Field;

    fn next(&mut self) -> Option<Field> {
        let tag = match self.reader.read_tag() {
            Ok(tag) => String::from_utf8(tag).unwrap(),
            Err(_) => return None,
        };

        let ty = self.reader.read_char().unwrap();

        let value = match ty {
            'A' => Value::Char(self.reader.read_char().unwrap()),
            'c' => Value::Int8(self.reader.read_i8().unwrap()),
            'C' => Value::UInt8(self.reader.read_u8().unwrap()),
            's' => Value::Int16(self.reader.read_i16().unwrap()),
            'S' => Value::UInt16(self.reader.read_u16().unwrap()),
            'i' => Value::Int32(self.reader.read_i32().unwrap()),
            'I' => Value::UInt32(self.reader.read_u32().unwrap()),
            'f' => Value::Float(self.reader.read_f32().unwrap()),
            'Z' => {
                let s = self.reader.read_string()
                    .map(|v| String::from_utf8(v).unwrap())
                    .unwrap();

                Value::String(s)
            },
            'H' => {
                let hex = self.reader.read_hex()
                    .map(|v| String::from_utf8(v).unwrap())
                    .unwrap();

                Value::Hex(hex)
            },
            'B' => {
                let ty = self.reader.read_char().unwrap();
                let len = self.reader.read_i32().unwrap() as usize;

                match ty {
                    'c' => Value::Int8Array(self.reader.read_i8_array(len).unwrap()),
                    'C' => Value::UInt8Array(self.reader.read_u8_array(len).unwrap()),
                    's' => Value::Int16Array(self.reader.read_i16_array(len).unwrap()),
                    'S' => Value::UInt16Array(self.reader.read_u16_array(len).unwrap()),
                    'i' => Value::Int32Array(self.reader.read_i32_array(len).unwrap()),
                    'I' => Value::UInt32Array(self.reader.read_u32_array(len).unwrap()),
                    'f' => Value::FloatArray(self.reader.read_f32_array(len).unwrap()),
                    _ => panic!("invalid field type '{}'", ty),
                }
            },
            _ => panic!("invalid field type '{}'", ty),
        };

        Some(Field::new(tag, value))
    }
}

#[derive(Debug)]
pub struct Data(Vec<u8>);

impl Data {
    pub fn new(raw_data: Vec<u8>) -> Data {
        Data(raw_data)
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn parse(&self) -> Vec<Field> {
        let cursor = Cursor::new(&self.0);
        let mut reader = Reader::new(cursor);
        reader.fields().collect()
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.0
    }

    pub fn extend(&mut self, field: &Field) {
        self.0.extend_from_slice(field.tag().as_bytes());

        self.0.push(field.value().ty() as u8);

        if let Some(subtype) = field.value().subtype() {
            self.0.push(subtype as u8);
        }

        match field.value() {
            Value::String(s) => {
                self.0.extend_from_slice(s.as_bytes());
                self.0.push(b'\0');
            },
            _ => unimplemented!(),
        }
    }
}

#[derive(Debug)]
pub enum Value {
    Char(char),
    Int8(i8),
    UInt8(u8),
    Int16(i16),
    UInt16(u16),
    Int32(i32),
    UInt32(u32),
    Float(f32),
    String(String),
    Hex(String),
    Int8Array(Vec<i8>),
    UInt8Array(Vec<u8>),
    Int16Array(Vec<i16>),
    UInt16Array(Vec<u16>),
    Int32Array(Vec<i32>),
    UInt32Array(Vec<u32>),
    FloatArray(Vec<f32>),
}

impl Value {
    pub fn ty(&self) -> char {
        match *self {
            Value::Char(_) => 'A',
            Value::Int8(_) => 'c',
            Value::UInt8(_) => 'C',
            Value::Int16(_) => 's',
            Value::UInt16(_) => 'S',
            Value::Int32(_) => 'i',
            Value::UInt32(_) => 'I',
            Value::Float(_) => 'f',
            Value::String(_) => 'Z',
            Value::Hex(_) => 'H',
            Value::Int8Array(_) => 'B',
            Value::UInt8Array(_) => 'B',
            Value::Int16Array(_) => 'B',
            Value::UInt16Array(_) => 'B',
            Value::Int32Array(_) => 'B',
            Value::UInt32Array(_) => 'B',
            Value::FloatArray(_) => 'B',
        }
    }

    pub fn subtype(&self) -> Option<char> {
        match *self {
            Value::Int8Array(_) => Some('c'),
            Value::UInt8Array(_) => Some('C'),
            Value::Int16Array(_) => Some('s'),
            Value::UInt16Array(_) => Some('S'),
            Value::Int32Array(_) => Some('i'),
            Value::UInt32Array(_) => Some('I'),
            Value::FloatArray(_) => Some('f'),
            _ => None,
        }
    }

    pub fn as_char(&self) -> Option<char> {
        match *self {
            Value::Char(c) => Some(c),
            _ => None,
        }
    }

    pub fn is_char(&self) -> bool {
        self.as_char().is_some()
    }

    pub fn as_str(&self) -> Option<&str> {
        match *self {
            Value::String(ref s) => Some(s),
            _ => None,
        }
    }

    pub fn is_str(&self) -> bool {
        self.as_str().is_some()
    }
}

#[derive(Debug)]
pub struct Field {
    tag: String,
    value: Value,
}

impl Field {
    pub fn new(tag: String, value: Value) -> Field {
        Field { tag, value }
    }

    pub fn tag(&self) -> &str {
        &self.tag
    }

    pub fn value(&self) -> &Value {
        &self.value
    }
}
