mod subtype;

use std::io::{self, Write};

use self::subtype::write_subtype;
use crate::{
    alignment::record::data::field::value::{array::Subtype, Array},
    writer::num,
};

pub(super) fn write_array<W>(writer: &mut W, array: &Array) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: u8 = b',';

    match array {
        Array::Int8(values) => {
            write_subtype(writer, Subtype::Int8)?;

            for result in values.iter() {
                let n = result?;
                writer.write_all(&[DELIMITER])?;
                num::write_i8(writer, n)?;
            }
        }
        Array::UInt8(values) => {
            write_subtype(writer, Subtype::UInt8)?;

            for result in values.iter() {
                let n = result?;
                writer.write_all(&[DELIMITER])?;
                num::write_u8(writer, n)?;
            }
        }
        Array::Int16(values) => {
            write_subtype(writer, Subtype::Int16)?;

            for result in values.iter() {
                let n = result?;
                writer.write_all(&[DELIMITER])?;
                num::write_i16(writer, n)?;
            }
        }
        Array::UInt16(values) => {
            write_subtype(writer, Subtype::UInt16)?;

            for result in values.iter() {
                let n = result?;
                writer.write_all(&[DELIMITER])?;
                num::write_u16(writer, n)?;
            }
        }
        Array::Int32(values) => {
            write_subtype(writer, Subtype::Int32)?;

            for result in values.iter() {
                let n = result?;
                writer.write_all(&[DELIMITER])?;
                num::write_i32(writer, n)?;
            }
        }
        Array::UInt32(values) => {
            write_subtype(writer, Subtype::UInt32)?;

            for result in values.iter() {
                let n = result?;
                writer.write_all(&[DELIMITER])?;
                num::write_u32(writer, n)?;
            }
        }
        Array::Float(values) => {
            write_subtype(writer, Subtype::Float)?;

            for result in values.iter() {
                let n = result?;
                write!(writer, ",{n}")?;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
pub(super) mod tests {
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

    impl<'a, N> crate::alignment::record::data::field::value::array::Values<'a, N> for T<'a, N>
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
            b"c,1,-2",
        )?;
        t(&mut buf, &Array::UInt8(Box::new(T::new(&[3, 5]))), b"C,3,5")?;
        t(
            &mut buf,
            &Array::Int16(Box::new(T::new(&[8, -13]))),
            b"s,8,-13",
        )?;
        t(
            &mut buf,
            &Array::UInt16(Box::new(T::new(&[21, 34]))),
            b"S,21,34",
        )?;
        t(
            &mut buf,
            &Array::Int32(Box::new(T::new(&[55, -89]))),
            b"i,55,-89",
        )?;
        t(
            &mut buf,
            &Array::UInt32(Box::new(T::new(&[144, 223]))),
            b"I,144,223",
        )?;
        t(
            &mut buf,
            &Array::Float(Box::new(T::new(&[8.0, 13.0]))),
            b"f,8,13",
        )?;

        Ok(())
    }
}
