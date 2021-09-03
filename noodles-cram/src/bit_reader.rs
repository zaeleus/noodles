use std::io::{self, Read};

use byteorder::ReadBytesExt;

pub struct BitReader<R>
where
    R: Read,
{
    inner: R,
    buf: u8,
    i: usize,
}

impl<R> BitReader<R>
where
    R: Read,
{
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            buf: 0,
            i: 8,
        }
    }

    pub fn read_u32(&mut self, n: usize) -> io::Result<u32> {
        let mut value = 0;

        for _ in 0..n {
            let bit = self.read_bit().map(u32::from)?;
            value = (value << 1) | bit;
        }

        Ok(value)
    }

    fn read_bit(&mut self) -> io::Result<u8> {
        if self.i >= 8 {
            self.buf = self.inner.read_u8()?;
            self.i = 0;
        }

        let bit = (self.buf >> (8 - self.i - 1)) & 0x01;

        self.i += 1;

        Ok(bit)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_u32() -> io::Result<()> {
        let data = [0b11001111, 0b01000000];
        let mut reader = BitReader::new(&data[..]);
        assert_eq!(reader.read_u32(4)?, 0x0c);
        assert_eq!(reader.read_u32(2)?, 0x03);
        assert_eq!(reader.read_u32(6)?, 0x34);
        Ok(())
    }
}
