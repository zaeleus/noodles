use std::io::{self, BufRead};

use byteorder::{LittleEndian, ReadBytesExt};

use super::{Field, Value};

pub struct Reader<R: BufRead> {
    reader: R,
}

impl<R: BufRead> Reader<R> {
    pub fn new(reader: R) -> Reader<R> {
        Reader { reader }
    }

    pub fn fields(&mut self) -> Fields<R> {
        Fields { reader: self }
    }

    fn read_field(&mut self) -> io::Result<Option<Field>> {
        let tag = match self.read_tag() {
            Ok(s) => {
                String::from_utf8(s).map_err(|e| {
                    io::Error::new(io::ErrorKind::InvalidInput, format!("{}", e))
                })?
            },
            Err(_) => return Ok(None),
        };

        let ty = self.read_char()?;

        let value = match ty {
            'A' => Value::Char(self.read_char()?),
            'c' => Value::Int8(self.read_i8()?),
            'C' => Value::UInt8(self.read_u8()?),
            's' => Value::Int16(self.read_i16()?),
            'S' => Value::UInt16(self.read_u16()?),
            'i' => Value::Int32(self.read_i32()?),
            'I' => Value::UInt32(self.read_u32()?),
            'f' => Value::Float(self.read_f32()?),
            'Z' => {
                self.read_string()
                    .and_then(|v| {
                        String::from_utf8(v).map_err(|e| {
                            io::Error::new(io::ErrorKind::InvalidInput, format!("{}", e))
                        })
                    })
                    .map(|v| Value::String(v))?
            },
            'H' => {
                self.read_hex()
                    .and_then(|v| {
                        String::from_utf8(v).map_err(|e| {
                            io::Error::new(io::ErrorKind::InvalidInput, format!("{}", e))
                        })
                    })
                    .map(|v| Value::Hex(v))?
            },
            'B' => {
                let ty = self.read_char()?;
                let len = self.read_i32()? as usize;

                match ty {
                    'c' => Value::Int8Array(self.read_i8_array(len)?),
                    'C' => Value::UInt8Array(self.read_u8_array(len)?),
                    's' => Value::Int16Array(self.read_i16_array(len)?),
                    'S' => Value::UInt16Array(self.read_u16_array(len)?),
                    'i' => Value::Int32Array(self.read_i32_array(len)?),
                    'I' => Value::UInt32Array(self.read_u32_array(len)?),
                    'f' => Value::FloatArray(self.read_f32_array(len)?),
                    _ => panic!("invalid field type 'B{}'", ty),
                }
            },
            _ => panic!("invalid field type '{}'", ty),
        };

        Ok(Some(Field::new(tag, value)))
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

pub struct Fields<'a, R: 'a + BufRead> {
    reader: &'a mut Reader<R>,
}

impl<'a, R: 'a + BufRead> Iterator for Fields<'a, R> {
    type Item = io::Result<Field>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_field() {
            Ok(Some(field)) => Some(Ok(field)),
            Ok(None) => None,
            Err(e) => Some(Err(e))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Reader;
    use crate::formats::bam::data::Field;

    #[test]
    fn test_read_field() {
        fn read_one_field(data: &[u8]) -> Field {
            let mut reader = Reader::new(&data[..]);
            reader.read_field().unwrap().unwrap()
        }

        let mut reader = Reader::new(&b""[..]);
        let field = reader.read_field().unwrap();
        assert!(field.is_none());

        let field = read_one_field(&b"WAAm"[..]);
        assert_eq!(field.tag(), "WA");
        assert!(field.value().is_char());
        assert_eq!(field.value().as_char(), Some('m'));

        let field = read_one_field(&b"JSc\xf9"[..]);
        assert_eq!(field.tag(), "JS");
        assert!(field.value().is_int8());
        assert_eq!(field.value().as_int8(), Some(-7));

        let field = read_one_field(&b"NSC\x07"[..]);
        assert_eq!(field.tag(), "NS");
        assert!(field.value().is_uint8());
        assert_eq!(field.value().as_uint8(), Some(7));

        let field = read_one_field(&b"UEs\xc0\xe0"[..]);
        assert_eq!(field.tag(), "UE");
        assert!(field.value().is_int16());
        assert_eq!(field.value().as_int16(), Some(-8000));

        let field = read_one_field(&b"QRS\x40\x1f"[..]);
        assert_eq!(field.tag(), "QR");
        assert!(field.value().is_uint16());
        assert_eq!(field.value().as_uint16(), Some(8000));

        let field = read_one_field(&b"DIi\x79\x96\xFF\xFF"[..]);
        assert_eq!(field.tag(), "DI");
        assert!(field.value().is_int32());
        assert_eq!(field.value().as_int32(), Some(-27015));

        let field = read_one_field(&b"TEI\x87\x69\x00\x00"[..]);
        assert_eq!(field.tag(), "TE");
        assert!(field.value().is_uint32());
        assert_eq!(field.value().as_uint32(), Some(27015));

        let field = read_one_field(&b"JBf\xd0\x0f\x49\x40"[..]);
        assert_eq!(field.tag(), "JB");
        assert!(field.value().is_float());
        let f = field.value().as_float().unwrap();
        assert!(f - 3.14159 < ::std::f32::EPSILON);

        let field = read_one_field(&b"NRZ\x6e\x6f\x6f\x64\x6c\x65\x73\x00"[..]);
        assert_eq!(field.tag(), "NR");
        assert!(field.value().is_str());
        assert_eq!(field.value().as_str(), Some("noodles"));

        let field = read_one_field(&b"LNH\x6e\x6f\x00"[..]);
        assert_eq!(field.tag(), "LN");
        assert!(field.value().is_hex());
        assert_eq!(field.value().as_hex(), Some("\x6e\x6f"));

        let field = read_one_field(&b"JHBc\x04\x00\x00\x00\xfd\xf9\xfc\xff"[..]);
        assert_eq!(field.tag(), "JH");
        assert!(field.value().is_int8_array());
        assert_eq!(field.value().as_int8_array(), Some(&[-3, -7, -4, -1][..]));

        let field = read_one_field(&b"LFBC\x04\x00\x00\x00\x03\x07\x04\x01"[..]);
        assert_eq!(field.tag(), "LF");
        assert!(field.value().is_uint8_array());
        assert_eq!(field.value().as_uint8_array(), Some(&[3, 7, 4, 1][..]));

        let field = read_one_field(&b"MZBs\x02\x00\x00\x00\xc0\xe0\xbf\xe0"[..]);
        assert_eq!(field.tag(), "MZ");
        assert!(field.value().is_int16_array());
        assert_eq!(field.value().as_int16_array(), Some(&[-8000, -8001][..]));

        let field = read_one_field(&b"JABS\x02\x00\x00\x00\x40\x1f\x41\x1f"[..]);
        assert_eq!(field.tag(), "JA");
        assert!(field.value().is_uint16_array());
        assert_eq!(field.value().as_uint16_array(), Some(&[8000, 8001][..]));

        let field = read_one_field(&b"RBBi\x01\x00\x00\x00\x79\x96\xff\xff"[..]);
        assert_eq!(field.tag(), "RB");
        assert!(field.value().is_int32_array());
        assert_eq!(field.value().as_int32_array(), Some(&[-27015][..]));

        let field = read_one_field(&b"IABI\x01\x00\x00\x00\x87\x69\x00\x00"[..]);
        assert_eq!(field.tag(), "IA");
        assert!(field.value().is_uint32_array());
        assert_eq!(field.value().as_uint32_array(), Some(&[27015][..]));

        let field = read_one_field(&b"EQBf\x01\x00\x00\x00\xa1\x67\x0b\x40"[..]);
        assert_eq!(field.tag(), "EQ");
        assert!(field.value().is_float_array());
        let a = field.value().as_float_array().unwrap();
        assert_eq!(a.len(), 1);
        assert!(a[0] - 2.7182 < ::std::f32::EPSILON);
    }
}
