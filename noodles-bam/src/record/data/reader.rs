//! BAM record data reader and iterators.

mod fields;

pub use self::fields::Fields;

use std::{
    convert::TryFrom,
    io::{self, BufRead},
};

use byteorder::{LittleEndian, ReadBytesExt};

use super::{
    field::{
        value::{Subtype, Type},
        Value,
    },
    Field,
};

/// A BAM record data reader.
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
    /// Creates a BAM record data reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::data::Reader;
    ///
    /// // NH:i:1  RG:Z:rg0
    /// let data = [
    ///     0x4e, 0x48, 0x69, 0x01, 0x00, 0x00, 0x00,
    ///     0x52, 0x47, 0x5a, 0x72, 0x67, 0x30, 0x00,
    /// ];
    /// let reader = Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a BAM record data field value.
    ///
    /// The stream is expected to be at the start of the value, i.e., after the tag and data type.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::data::{field::{value::Type, Value}, Reader};
    ///
    /// let data = [0x01, 0x00, 0x00, 0x00];
    /// let mut reader = Reader::new(&data[..]);
    ///
    /// let value = reader.read_value_type(Type::Int32)?;
    ///
    /// assert_eq!(value, Value::Int32(1));
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_value_type(&mut self, ty: Type) -> io::Result<Value> {
        read_value_type(&mut self.inner, ty)
    }

    /// Returns an iterator over data fields.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::{data::{field::Value, Field, Reader}};
    /// use noodles_sam::record::data::field::Tag;
    ///
    /// // NH:i:1  RG:Z:rg0
    /// let data = [
    ///     0x4e, 0x48, 0x69, 0x01, 0x00, 0x00, 0x00,
    ///     0x52, 0x47, 0x5a, 0x72, 0x67, 0x30, 0x00,
    /// ];
    /// let reader = Reader::new(&data[..]);
    ///
    /// let mut fields = reader.fields();
    ///
    /// let field = fields.next().transpose()?;
    /// assert_eq!(field, Some(Field::new(Tag::AlignmentHitCount, Value::Int32(1))));
    ///
    /// let field = fields.next().transpose()?;
    /// assert_eq!(field, Some(Field::new(Tag::ReadGroup, Value::String(String::from("rg0")))));
    ///
    /// assert!(fields.next().is_none());
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn fields(self) -> Fields<R> {
        Fields::new(self)
    }

    fn read_field(&mut self) -> io::Result<Option<Field>> {
        use crate::reader::data::field::{tag::read_tag, value::ty::read_type};

        let tag = match read_tag(&mut self.inner) {
            Ok(t) => t,
            Err(_) => return Ok(None),
        };

        let ty = read_type(&mut self.inner)?;
        let value = read_value_type(&mut self.inner, ty)?;

        Ok(Some(Field::new(tag, value)))
    }
}

fn read_value_type<R>(reader: &mut R, ty: Type) -> io::Result<Value>
where
    R: BufRead,
{
    match ty {
        Type::Char => reader.read_u8().map(char::from).map(Value::Char),
        Type::Int8 => reader.read_i8().map(Value::Int8),
        Type::UInt8 => reader.read_u8().map(Value::UInt8),
        Type::Int16 => reader.read_i16::<LittleEndian>().map(Value::Int16),
        Type::UInt16 => reader.read_u16::<LittleEndian>().map(Value::UInt16),
        Type::Int32 => reader.read_i32::<LittleEndian>().map(Value::Int32),
        Type::UInt32 => reader.read_u32::<LittleEndian>().map(Value::UInt32),
        Type::Float => reader.read_f32::<LittleEndian>().map(Value::Float),
        Type::String => read_string(reader).map(Value::String),
        Type::Hex => read_string(reader).map(Value::Hex),
        Type::Array => read_array(reader),
    }
}

fn read_string<R>(reader: &mut R) -> io::Result<String>
where
    R: BufRead,
{
    let mut buf = Vec::new();
    reader.read_until(b'\0', &mut buf)?;
    buf.pop();
    String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_array<R>(reader: &mut R) -> io::Result<Value>
where
    R: BufRead,
{
    let subtype = reader.read_u8().and_then(|b| {
        Subtype::try_from(b).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let len = reader.read_i32::<LittleEndian>()? as usize;

    match subtype {
        Subtype::Int8 => {
            let mut buf = vec![0; len];
            reader.read_i8_into(&mut buf)?;
            Ok(Value::Int8Array(buf))
        }
        Subtype::UInt8 => {
            let mut buf = vec![0; len];
            reader.read_exact(&mut buf)?;
            Ok(Value::UInt8Array(buf))
        }
        Subtype::Int16 => {
            let mut buf = vec![0; len];
            reader.read_i16_into::<LittleEndian>(&mut buf)?;
            Ok(Value::Int16Array(buf))
        }
        Subtype::UInt16 => {
            let mut buf = vec![0; len];
            reader.read_u16_into::<LittleEndian>(&mut buf)?;
            Ok(Value::UInt16Array(buf))
        }
        Subtype::Int32 => {
            let mut buf = vec![0; len];
            reader.read_i32_into::<LittleEndian>(&mut buf)?;
            Ok(Value::Int32Array(buf))
        }
        Subtype::UInt32 => {
            let mut buf = vec![0; len];
            reader.read_u32_into::<LittleEndian>(&mut buf)?;
            Ok(Value::UInt32Array(buf))
        }
        Subtype::Float => {
            let mut buf = vec![0.0; len];
            reader.read_f32_into::<LittleEndian>(&mut buf)?;
            Ok(Value::FloatArray(buf))
        }
    }
}

#[cfg(test)]
mod tests {
    use noodles_sam::record::data::field::Tag;

    use super::*;

    #[test]
    fn test_read_field() {
        fn read_one_field(data: &[u8]) -> Field {
            let mut reader = Reader::new(data);
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
