mod subtype;

use std::io;

use noodles_sam::alignment::record::data::field::value::{Array, array::Subtype};

use crate::record::codec::encoder::num::{
    write_f32_le, write_i8, write_i16_le, write_i32_le, write_u8, write_u16_le, write_u32_le,
};

use self::subtype::write_subtype;

pub fn write_array(dst: &mut Vec<u8>, array: &Array) -> io::Result<()> {
    match array {
        Array::Int8(values) => {
            write_header(dst, Subtype::Int8, values.len())?;

            for result in values.iter() {
                let n = result?;
                write_i8(dst, n);
            }
        }
        Array::UInt8(values) => {
            write_header(dst, Subtype::UInt8, values.len())?;

            for result in values.iter() {
                let n = result?;
                write_u8(dst, n);
            }
        }
        Array::Int16(values) => {
            write_header(dst, Subtype::Int16, values.len())?;

            for result in values.iter() {
                let n = result?;
                write_i16_le(dst, n);
            }
        }
        Array::UInt16(values) => {
            write_header(dst, Subtype::UInt16, values.len())?;

            for result in values.iter() {
                let n = result?;
                write_u16_le(dst, n);
            }
        }
        Array::Int32(values) => {
            write_header(dst, Subtype::Int32, values.len())?;

            for result in values.iter() {
                let n = result?;
                write_i32_le(dst, n);
            }
        }
        Array::UInt32(values) => {
            write_header(dst, Subtype::UInt32, values.len())?;

            for result in values.iter() {
                let n = result?;
                write_u32_le(dst, n);
            }
        }
        Array::Float(values) => {
            write_header(dst, Subtype::Float, values.len())?;

            for result in values.iter() {
                let n = result?;
                write_f32_le(dst, n);
            }
        }
    }

    Ok(())
}

pub fn write_header(dst: &mut Vec<u8>, subtype: Subtype, len: usize) -> io::Result<()> {
    write_subtype(dst, subtype);

    let n = u32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u32_le(dst, n);

    Ok(())
}

#[cfg(test)]
pub(super) mod tests {
    use noodles_sam as sam;

    use super::*;

    pub(crate) struct T<'a, N>(&'a [N]);

    impl<'a, N> T<'a, N>
    where
        N: Copy,
    {
        pub(crate) fn new(values: &'a [N]) -> Self {
            Self(values)
        }
    }

    impl<'a, N> sam::alignment::record::data::field::value::array::Values<'a, N> for T<'a, N>
    where
        N: Copy,
    {
        fn len(&self) -> usize {
            self.0.len()
        }

        fn iter(&self) -> Box<dyn Iterator<Item = io::Result<N>> + '_> {
            Box::new(self.0.iter().copied().map(Ok))
        }
    }

    #[test]
    fn test_write_array() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, array: &Array, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_array(buf, array)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(
            &mut buf,
            &Array::Int8(Box::new(T::new(&[1, -2]))),
            &[
                b'c', // subtype = Int8
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x01, // values[0] = 1
                0xfe, // values[1] = -2
            ],
        )?;

        t(
            &mut buf,
            &Array::UInt8(Box::new(T::new(&[3, 5]))),
            &[
                b'C', // subtype = UInt8
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x03, // values[0] = 3
                0x05, // values[1] = 5
            ],
        )?;

        t(
            &mut buf,
            &Array::Int16(Box::new(T::new(&[8, -13]))),
            &[
                b's', // subtype = Int16
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x08, 0x00, // values[0] = 8
                0xf3, 0xff, // values[1] = -13
            ],
        )?;

        t(
            &mut buf,
            &Array::UInt16(Box::new(T::new(&[21, 34]))),
            &[
                b'S', // subtype = UInt16
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x15, 0x00, // values[0] = 21
                0x22, 0x00, // values[1] = 34
            ],
        )?;

        t(
            &mut buf,
            &Array::Int32(Box::new(T::new(&[55, -89]))),
            &[
                b'i', // subtype = Int32
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x37, 0x00, 0x00, 0x00, // values[0] = 55
                0xa7, 0xff, 0xff, 0xff, // values[1] = -89
            ],
        )?;

        t(
            &mut buf,
            &Array::UInt32(Box::new(T::new(&[144, 223]))),
            &[
                b'I', // subtype = UInt32
                0x02, 0x00, 0x00, 0x00, // count = 2
                0x90, 0x00, 0x00, 0x00, // values[0] = 55
                0xdf, 0x00, 0x00, 0x00, // values[1] = -89
            ],
        )?;

        t(
            &mut buf,
            &Array::Float(Box::new(T::new(&[8.0, 13.0]))),
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
