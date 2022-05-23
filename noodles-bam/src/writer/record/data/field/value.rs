mod subtype;
mod ty;

pub use self::{subtype::put_subtype, ty::put_type};

use std::{ffi::CString, io};

use bytes::BufMut;
use noodles_sam::record::data::field::{value::Subtype, Value};

/// Writes a BAM record data field value.
///
/// # Examples
///
/// ```
/// # use std::io;
/// use noodles_bam::writer::record::data::field::put_value;
/// use noodles_sam::record::data::field::Value;
/// let mut buf = Vec::new();
/// let value = Value::UInt8(0);
/// put_value(&mut buf, &value)?;
/// assert_eq!(buf, [0x00]);
/// # Ok::<_, io::Error>(())
/// ```
pub fn put_value<B>(dst: &mut B, value: &Value) -> io::Result<()>
where
    B: BufMut,
{
    match value {
        Value::Character(c) => dst.put_u8(u8::from(*c)),
        Value::Int8(n) => dst.put_i8(*n),
        Value::UInt8(n) => dst.put_u8(*n),
        Value::Int16(n) => dst.put_i16_le(*n),
        Value::UInt16(n) => dst.put_u16_le(*n),
        Value::Int32(n) => dst.put_i32_le(*n),
        Value::UInt32(n) => dst.put_u32_le(*n),
        Value::Float(n) => dst.put_f32_le(*n),
        Value::String(s) | Value::Hex(s) => put_string(dst, s)?,
        Value::Int8Array(values) => {
            put_array_header(dst, Subtype::Int8, values.len())?;

            for &n in values {
                dst.put_i8(n);
            }
        }
        Value::UInt8Array(values) => {
            put_array_header(dst, Subtype::UInt8, values.len())?;

            for &n in values {
                dst.put_u8(n);
            }
        }
        Value::Int16Array(values) => {
            put_array_header(dst, Subtype::Int16, values.len())?;

            for &n in values {
                dst.put_i16_le(n);
            }
        }
        Value::UInt16Array(values) => {
            put_array_header(dst, Subtype::UInt16, values.len())?;

            for &n in values {
                dst.put_u16_le(n);
            }
        }
        Value::Int32Array(values) => {
            put_array_header(dst, Subtype::Int32, values.len())?;

            for &n in values {
                dst.put_i32_le(n);
            }
        }
        Value::UInt32Array(values) => {
            put_array_header(dst, Subtype::UInt32, values.len())?;

            for &n in values {
                dst.put_u32_le(n);
            }
        }
        Value::FloatArray(values) => {
            put_array_header(dst, Subtype::Float, values.len())?;

            for &n in values {
                dst.put_f32_le(n);
            }
        }
    }

    Ok(())
}

fn put_string<B>(dst: &mut B, s: &str) -> io::Result<()>
where
    B: BufMut,
{
    let c_str = CString::new(s).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.put(c_str.as_bytes_with_nul());
    Ok(())
}

fn put_array_header<B>(dst: &mut B, subtype: Subtype, len: usize) -> io::Result<()>
where
    B: BufMut,
{
    put_subtype(dst, subtype);

    let n = u32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.put_u32_le(n);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_value() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_sam::record::data::field::value::Character;

        fn t(buf: &mut Vec<u8>, value: &Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            put_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(
            &mut buf,
            &Value::Character(Character::try_from('n')?),
            &[b'n'],
        )?;
        t(&mut buf, &Value::Int8(1), &[0x01])?;
        t(&mut buf, &Value::UInt8(2), &[0x02])?;
        t(&mut buf, &Value::Int16(3), &[0x03, 0x00])?;
        t(&mut buf, &Value::UInt16(5), &[0x05, 0x00])?;
        t(&mut buf, &Value::Int32(8), &[0x08, 0x00, 0x00, 0x00])?;
        t(&mut buf, &Value::UInt32(13), &[0x0d, 0x00, 0x00, 0x00])?;
        t(&mut buf, &Value::Float(8.0), &[0x00, 0x00, 0x00, 0x41])?;

        t(
            &mut buf,
            &Value::String(String::from("ndls")),
            &[b'n', b'd', b'l', b's', 0x00],
        )?;

        t(
            &mut buf,
            &Value::Hex(String::from("CAFE")),
            &[b'C', b'A', b'F', b'E', 0x00],
        )?;

        t(
            &mut buf,
            &Value::Int8Array(vec![1, -2]),
            &[
                b'c', // subtype = Int8
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x01, // values[0] = 1
                0xfe, // values[1] = -2
            ],
        )?;

        t(
            &mut buf,
            &Value::UInt8Array(vec![3, 5]),
            &[
                b'C', // subtype = UInt8
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x03, // values[0] = 3
                0x05, // values[1] = 5
            ],
        )?;

        t(
            &mut buf,
            &Value::Int16Array(vec![8, -13]),
            &[
                b's', // subtype = Int16
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x08, 0x00, // values[0] = 8
                0xf3, 0xff, // values[1] = -13
            ],
        )?;

        t(
            &mut buf,
            &Value::UInt16Array(vec![21, 34]),
            &[
                b'S', // subtype = UInt16
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x15, 0x00, // values[0] = 21
                0x22, 0x00, // values[1] = 34
            ],
        )?;

        t(
            &mut buf,
            &Value::Int32Array(vec![55, -89]),
            &[
                b'i', // subtype = Int32
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x37, 0x00, 0x00, 0x00, // values[0] = 55
                0xa7, 0xff, 0xff, 0xff, // values[1] = -89
            ],
        )?;

        t(
            &mut buf,
            &Value::UInt32Array(vec![144, 223]),
            &[
                b'I', // subtype = UInt32
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x90, 0x00, 0x00, 0x00, // values[0] = 55
                0xdf, 0x00, 0x00, 0x00, // values[1] = -89
            ],
        )?;

        t(
            &mut buf,
            &Value::FloatArray(vec![8.0, 13.0]),
            &[
                b'f', // subtype = Float
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x00, 0x00, 0x00, 0x41, // values[0] = 8.0
                0x00, 0x00, 0x50, 0x41, // values[1] = 13.0
            ],
        )?;

        Ok(())
    }
}
