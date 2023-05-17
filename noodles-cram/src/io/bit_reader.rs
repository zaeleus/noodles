use std::io;

use bytes::Buf;

pub struct BitReader<B>
where
    B: Buf,
{
    src: B,
    buf: u8,
    i: usize,
}

impl<B> BitReader<B>
where
    B: Buf,
{
    pub fn new(src: B) -> Self {
        Self { src, buf: 0, i: 8 }
    }

    pub fn read_u32(&mut self, n: u32) -> io::Result<u32> {
        let mut value = 0;

        for _ in 0..n {
            let bit = self.read_bit().map(u32::from)?;
            value = (value << 1) | bit;
        }

        Ok(value)
    }

    pub(crate) fn read_bit(&mut self) -> io::Result<u8> {
        if self.i >= 8 {
            if !self.src.has_remaining() {
                return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
            }

            self.buf = self.src.get_u8();
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
        assert_eq!(reader.read_u32(4)?, 0b1100);
        assert_eq!(reader.read_u32(2)?, 0b11);
        assert_eq!(reader.read_u32(6)?, 0b110100);
        Ok(())
    }
}
