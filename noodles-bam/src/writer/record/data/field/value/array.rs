use std::io;

use bytes::BufMut;
use noodles_sam::record::data::field::value::{Array, Subtype};

use super::put_subtype;

pub fn put_array<B>(dst: &mut B, array: &Array) -> io::Result<()>
where
    B: BufMut,
{
    match array {
        Array::Int8(values) => {
            put_header(dst, Subtype::Int8, values.len())?;

            for &n in values {
                dst.put_i8(n);
            }
        }
        Array::UInt8(values) => {
            put_header(dst, Subtype::UInt8, values.len())?;

            for &n in values {
                dst.put_u8(n);
            }
        }
        Array::Int16(values) => {
            put_header(dst, Subtype::Int16, values.len())?;

            for &n in values {
                dst.put_i16_le(n);
            }
        }
        Array::UInt16(values) => {
            put_header(dst, Subtype::UInt16, values.len())?;

            for &n in values {
                dst.put_u16_le(n);
            }
        }
        Array::Int32(values) => {
            put_header(dst, Subtype::Int32, values.len())?;

            for &n in values {
                dst.put_i32_le(n);
            }
        }
        Array::UInt32(values) => {
            put_header(dst, Subtype::UInt32, values.len())?;

            for &n in values {
                dst.put_u32_le(n);
            }
        }
        Array::Float(values) => {
            put_header(dst, Subtype::Float, values.len())?;

            for &n in values {
                dst.put_f32_le(n);
            }
        }
    }

    Ok(())
}

pub fn put_header<B>(dst: &mut B, subtype: Subtype, len: usize) -> io::Result<()>
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
    fn test_put_array() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, array: &Array, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            put_array(buf, array)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(
            &mut buf,
            &Array::Int8(vec![1, -2]),
            &[
                b'c', // subtype = Int8
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x01, // values[0] = 1
                0xfe, // values[1] = -2
            ],
        )?;

        t(
            &mut buf,
            &Array::UInt8(vec![3, 5]),
            &[
                b'C', // subtype = UInt8
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x03, // values[0] = 3
                0x05, // values[1] = 5
            ],
        )?;

        t(
            &mut buf,
            &Array::Int16(vec![8, -13]),
            &[
                b's', // subtype = Int16
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x08, 0x00, // values[0] = 8
                0xf3, 0xff, // values[1] = -13
            ],
        )?;

        t(
            &mut buf,
            &Array::UInt16(vec![21, 34]),
            &[
                b'S', // subtype = UInt16
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x15, 0x00, // values[0] = 21
                0x22, 0x00, // values[1] = 34
            ],
        )?;

        t(
            &mut buf,
            &Array::Int32(vec![55, -89]),
            &[
                b'i', // subtype = Int32
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x37, 0x00, 0x00, 0x00, // values[0] = 55
                0xa7, 0xff, 0xff, 0xff, // values[1] = -89
            ],
        )?;

        t(
            &mut buf,
            &Array::UInt32(vec![144, 223]),
            &[
                b'I', // subtype = UInt32
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x90, 0x00, 0x00, 0x00, // values[0] = 55
                0xdf, 0x00, 0x00, 0x00, // values[1] = -89
            ],
        )?;

        t(
            &mut buf,
            &Array::Float(vec![8.0, 13.0]),
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
