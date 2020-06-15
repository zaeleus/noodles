mod fields;

pub use self::fields::Fields;

use std::io::{self, BufRead};

use byteorder::{LittleEndian, ReadBytesExt};

use super::{field::Value, Field};

pub struct Reader<R>
where
    R: BufRead,
{
    inner: R,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    pub fn fields(self) -> Fields<R> {
        Fields::new(self)
    }

    fn read_field(&mut self) -> io::Result<Option<Field>> {
        let tag = match self.read_tag() {
            Ok(data) => String::from_utf8(data)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                .and_then(|s| {
                    s.parse()
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                })?,
            Err(_) => return Ok(None),
        };

        let value = self.read_value()?;

        Ok(Some(Field::new(tag, value)))
    }

    fn read_tag(&mut self) -> io::Result<Vec<u8>> {
        let mut buf = vec![0; 2];
        self.inner.read_exact(&mut buf)?;
        Ok(buf)
    }

    fn read_value(&mut self) -> io::Result<Value> {
        let ty = self.read_char()?;

        match ty {
            'A' => self.read_char().map(Value::Char),
            'c' => self.read_i8().map(Value::Int8),
            'C' => self.read_u8().map(Value::UInt8),
            's' => self.read_i16().map(Value::Int16),
            'S' => self.read_u16().map(Value::UInt16),
            'i' => self.read_i32().map(Value::Int32),
            'I' => self.read_u32().map(Value::UInt32),
            'f' => self.read_f32().map(Value::Float),
            'Z' => self
                .read_string()
                .and_then(|v| {
                    String::from_utf8(v).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                })
                .map(Value::String),
            'H' => self
                .read_hex()
                .and_then(|v| {
                    String::from_utf8(v).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                })
                .map(Value::Hex),
            'B' => {
                let ty = self.read_char()?;
                let len = self.read_i32()? as usize;

                match ty {
                    'c' => self.read_i8_array(len).map(Value::Int8Array),
                    'C' => self.read_u8_array(len).map(Value::UInt8Array),
                    's' => self.read_i16_array(len).map(Value::Int16Array),
                    'S' => self.read_u16_array(len).map(Value::UInt16Array),
                    'i' => self.read_i32_array(len).map(Value::Int32Array),
                    'I' => self.read_u32_array(len).map(Value::UInt32Array),
                    'f' => self.read_f32_array(len).map(Value::FloatArray),
                    _ => Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid field type 'B{}'", ty),
                    )),
                }
            }
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid field type '{}'", ty),
            )),
        }
    }

    fn read_char(&mut self) -> io::Result<char> {
        self.inner.read_u8().map(|b| b as char)
    }

    fn read_i8(&mut self) -> io::Result<i8> {
        self.inner.read_i8()
    }

    fn read_u8(&mut self) -> io::Result<u8> {
        self.inner.read_u8()
    }

    fn read_i16(&mut self) -> io::Result<i16> {
        self.inner.read_i16::<LittleEndian>()
    }

    fn read_u16(&mut self) -> io::Result<u16> {
        self.inner.read_u16::<LittleEndian>()
    }

    fn read_i32(&mut self) -> io::Result<i32> {
        self.inner.read_i32::<LittleEndian>()
    }

    fn read_u32(&mut self) -> io::Result<u32> {
        self.inner.read_u32::<LittleEndian>()
    }

    fn read_f32(&mut self) -> io::Result<f32> {
        self.inner.read_f32::<LittleEndian>()
    }

    fn read_string(&mut self) -> io::Result<Vec<u8>> {
        let mut buf = Vec::new();
        self.inner.read_until(b'\0', &mut buf)?;
        buf.pop();
        Ok(buf)
    }

    fn read_hex(&mut self) -> io::Result<Vec<u8>> {
        self.read_string()
    }

    fn read_i8_array(&mut self, len: usize) -> io::Result<Vec<i8>> {
        let mut buf = vec![0; len];
        self.inner.read_exact(&mut buf)?;
        Ok(buf.iter().map(|&b| b as i8).collect())
    }

    fn read_u8_array(&mut self, len: usize) -> io::Result<Vec<u8>> {
        let mut buf = vec![0; len];
        self.inner.read_exact(&mut buf)?;
        Ok(buf)
    }

    fn read_i16_array(&mut self, len: usize) -> io::Result<Vec<i16>> {
        let mut buf = vec![0; len];
        self.inner.read_i16_into::<LittleEndian>(&mut buf)?;
        Ok(buf)
    }

    fn read_u16_array(&mut self, len: usize) -> io::Result<Vec<u16>> {
        let mut buf = vec![0; len];
        self.inner.read_u16_into::<LittleEndian>(&mut buf)?;
        Ok(buf)
    }

    fn read_i32_array(&mut self, len: usize) -> io::Result<Vec<i32>> {
        let mut buf = vec![0; len];
        self.inner.read_i32_into::<LittleEndian>(&mut buf)?;
        Ok(buf)
    }

    fn read_u32_array(&mut self, len: usize) -> io::Result<Vec<u32>> {
        let mut buf = vec![0; len];
        self.inner.read_u32_into::<LittleEndian>(&mut buf)?;
        Ok(buf)
    }

    fn read_f32_array(&mut self, len: usize) -> io::Result<Vec<f32>> {
        let mut buf = vec![0.0; len];
        self.inner.read_f32_into::<LittleEndian>(&mut buf)?;
        Ok(buf)
    }
}

#[cfg(test)]
mod tests {
    use noodles_sam::record::data::field::Tag;

    use super::*;

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
        assert_eq!(field.tag(), &Tag::Other(String::from("WA")));
        assert!(field.value().is_char());
        assert_eq!(field.value().as_char(), Some('m'));

        let field = read_one_field(&b"JSc\xf9"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("JS")));
        assert!(field.value().is_int8());
        assert_eq!(field.value().as_int8(), Some(-7));

        let field = read_one_field(&b"NSC\x07"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("NS")));
        assert!(field.value().is_uint8());
        assert_eq!(field.value().as_uint8(), Some(7));

        let field = read_one_field(&b"UEs\xc0\xe0"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("UE")));
        assert!(field.value().is_int16());
        assert_eq!(field.value().as_int16(), Some(-8000));

        let field = read_one_field(&b"QRS\x40\x1f"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("QR")));
        assert!(field.value().is_uint16());
        assert_eq!(field.value().as_uint16(), Some(8000));

        let field = read_one_field(&b"DIi\x79\x96\xFF\xFF"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("DI")));
        assert!(field.value().is_int32());
        assert_eq!(field.value().as_int32(), Some(-27015));

        let field = read_one_field(&b"TEI\x87\x69\x00\x00"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("TE")));
        assert!(field.value().is_uint32());
        assert_eq!(field.value().as_uint32(), Some(27015));

        let field = read_one_field(&b"JBf\x56\x0e\x49\x40"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("JB")));
        assert!(field.value().is_float());
        let f = field.value().as_float().unwrap();
        assert_eq!(f, 3.1415);

        let field = read_one_field(&b"NRZ\x6e\x6f\x6f\x64\x6c\x65\x73\x00"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("NR")));
        assert!(field.value().is_str());
        assert_eq!(field.value().as_str(), Some("noodles"));

        let field = read_one_field(&b"LNH\x6e\x6f\x00"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("LN")));
        assert!(field.value().is_hex());
        assert_eq!(field.value().as_hex(), Some("\x6e\x6f"));

        let field = read_one_field(&b"JHBc\x04\x00\x00\x00\xfd\xf9\xfc\xff"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("JH")));
        assert!(field.value().is_int8_array());
        assert_eq!(field.value().as_int8_array(), Some(&[-3, -7, -4, -1][..]));

        let field = read_one_field(&b"LFBC\x04\x00\x00\x00\x03\x07\x04\x01"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("LF")));
        assert!(field.value().is_uint8_array());
        assert_eq!(field.value().as_uint8_array(), Some(&[3, 7, 4, 1][..]));

        let field = read_one_field(&b"MZBs\x02\x00\x00\x00\xc0\xe0\xbf\xe0"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("MZ")));
        assert!(field.value().is_int16_array());
        assert_eq!(field.value().as_int16_array(), Some(&[-8000, -8001][..]));

        let field = read_one_field(&b"JABS\x02\x00\x00\x00\x40\x1f\x41\x1f"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("JA")));
        assert!(field.value().is_uint16_array());
        assert_eq!(field.value().as_uint16_array(), Some(&[8000, 8001][..]));

        let field = read_one_field(&b"RBBi\x01\x00\x00\x00\x79\x96\xff\xff"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("RB")));
        assert!(field.value().is_int32_array());
        assert_eq!(field.value().as_int32_array(), Some(&[-27015][..]));

        let field = read_one_field(&b"IABI\x01\x00\x00\x00\x87\x69\x00\x00"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("IA")));
        assert!(field.value().is_uint32_array());
        assert_eq!(field.value().as_uint32_array(), Some(&[27015][..]));

        let field = read_one_field(&b"EQBf\x01\x00\x00\x00\xa1\xf8\x2d\x40"[..]);
        assert_eq!(field.tag(), &Tag::Other(String::from("EQ")));
        assert!(field.value().is_float_array());
        let a = field.value().as_float_array().unwrap();
        assert_eq!(a.len(), 1);
        assert_eq!(a[0], 2.7183);
    }
}
