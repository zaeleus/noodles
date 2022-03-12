//! BAM record data field and values.

use std::{io, num};

use bytes::BufMut;
use noodles_sam::{
    self as sam,
    record::data::field::{Tag, Value},
};

/// A BAM record data field.
#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    tag: Tag,
    value: Value,
}

impl Field {
    /// Creates a data field.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::data::Field;
    /// use noodles_sam::record::data::field::{Tag, Value};
    /// let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// ```
    pub fn new(tag: Tag, value: Value) -> Self {
        Self { tag, value }
    }

    /// Returns the data field tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::data::Field;
    /// use noodles_sam::record::data::field::{Tag, Value};
    /// let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// assert_eq!(field.tag(), Tag::AlignmentHitCount);
    /// ```
    pub fn tag(&self) -> Tag {
        self.tag
    }

    /// Returns the data field value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::data::Field;
    /// use noodles_sam::record::data::field::{Tag, Value};
    /// let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));
    /// assert_eq!(field.value(), &Value::Int32(1));
    /// ```
    pub fn value(&self) -> &Value {
        &self.value
    }
}

impl From<Field> for sam::record::data::Field {
    fn from(field: Field) -> Self {
        Self::new(field.tag, field.value)
    }
}

impl TryFrom<Field> for Vec<u8> {
    type Error = io::Error;

    fn try_from(field: Field) -> Result<Self, Self::Error> {
        let mut buf = Vec::new();
        put_field(&mut buf, &field)?;
        Ok(buf)
    }
}

fn put_field<B>(buf: &mut B, field: &Field) -> io::Result<()>
where
    B: BufMut,
{
    put_field_tag(buf, field.tag());
    put_field_value_type(buf, field.value());
    put_field_value(buf, field.value())?;
    Ok(())
}

fn put_field_tag<B>(buf: &mut B, tag: Tag)
where
    B: BufMut,
{
    buf.put(&tag.as_ref()[..]);
}

fn put_field_value_type<B>(buf: &mut B, value: &Value)
where
    B: BufMut,
{
    buf.put_u8(u8::from(value.ty()));

    if let Some(subtype) = value.subtype() {
        buf.put_u8(u8::from(subtype));
    }
}

fn put_field_value<B>(buf: &mut B, value: &Value) -> io::Result<()>
where
    B: BufMut,
{
    fn invalid_array_len(e: num::TryFromIntError) -> io::Error {
        io::Error::new(io::ErrorKind::InvalidInput, e)
    }

    match value {
        Value::Char(c) => buf.put_u8(*c as u8),
        Value::Int8(n) => buf.put_i8(*n),
        Value::UInt8(n) => buf.put_u8(*n),
        Value::Int16(n) => buf.put_i16_le(*n),
        Value::UInt16(n) => buf.put_u16_le(*n),
        Value::Int32(n) => buf.put_i32_le(*n),
        Value::UInt32(n) => buf.put_u32_le(*n),
        Value::Float(n) => buf.put_f32_le(*n),
        Value::String(s) | Value::Hex(s) => put_field_string_value(buf, s)?,
        Value::Int8Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            buf.put_u32_le(len);

            for &n in values {
                buf.put_i8(n);
            }
        }
        Value::UInt8Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            buf.put_u32_le(len);

            for &n in values {
                buf.put_u8(n);
            }
        }
        Value::Int16Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            buf.put_u32_le(len);

            for &n in values {
                buf.put_i16_le(n);
            }
        }
        Value::UInt16Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            buf.put_u32_le(len);

            for &n in values {
                buf.put_u16_le(n);
            }
        }
        Value::Int32Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            buf.put_u32_le(len);

            for &n in values {
                buf.put_i32_le(n);
            }
        }
        Value::UInt32Array(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            buf.put_u32_le(len);

            for &n in values {
                buf.put_u32_le(n);
            }
        }
        Value::FloatArray(values) => {
            let len = u32::try_from(values.len()).map_err(invalid_array_len)?;
            buf.put_u32_le(len);

            for &n in values {
                buf.put_f32_le(n);
            }
        }
    }

    Ok(())
}

fn put_field_string_value<B>(buf: &mut B, s: &str) -> io::Result<()>
where
    B: BufMut,
{
    use std::ffi::CString;

    let c_str = CString::new(s).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    buf.put(c_str.as_bytes_with_nul());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_field_for_sam_record_data_field() {
        let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));

        let actual = sam::record::data::Field::from(field);
        let expected = sam::record::data::Field::new(
            Tag::AlignmentHitCount,
            sam::record::data::field::Value::Int32(1),
        );

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_try_from_field_for_vec_u8() -> io::Result<()> {
        let field = Field::new(Tag::AlignmentHitCount, Value::Int32(1));

        let actual = <Vec<u8>>::try_from(field)?;
        let expected = vec![b'N', b'H', b'i', 0x01, 0x00, 0x00, 0x00];

        assert_eq!(actual, expected);

        Ok(())
    }
}
