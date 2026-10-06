use std::io::{self, Write};

use lexical_core::FormattedSize;

pub fn write_i8<W>(writer: &mut W, n: i8) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; i8::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_u8<W>(writer: &mut W, n: u8) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; u8::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_i16<W>(writer: &mut W, n: i16) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; i16::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_u16<W>(writer: &mut W, n: u16) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; u16::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_i32<W>(writer: &mut W, n: i32) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; i32::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_u32<W>(writer: &mut W, n: u32) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; u32::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_usize<W>(writer: &mut W, n: usize) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; usize::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_f32<W>(writer: &mut W, n: f32) -> io::Result<()>
where
    W: Write,
{
    const FORMAT: u128 = lexical_core::format::STANDARD;

    const OPTIONS: lexical_core::WriteFloatOptions = lexical_core::WriteFloatOptionsBuilder::new()
        .trim_floats(true)
        .build_strict();

    let mut dst = [0; f32::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write_with_options::<_, FORMAT>(n, &mut dst, &OPTIONS);
    writer.write_all(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_f32() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, n: f32, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_f32(buf, n)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, -1.0, b"-1")?;
        t(&mut buf, -0.5, b"-0.5")?;
        t(&mut buf, 0.0, b"0")?;
        t(&mut buf, 0.5, b"0.5")?;
        t(&mut buf, 1.0, b"1")?;

        Ok(())
    }
}
