use std::io;

pub struct BitReader<'c> {
    src: &'c [u8],
    buf: u8,
    i: usize,
}

impl<'c> BitReader<'c> {
    pub fn new(src: &'c [u8]) -> Self {
        Self {
            src,
            buf: 0x00,
            i: 8,
        }
    }

    pub fn read_u32(&mut self, len: u32) -> io::Result<u32> {
        let mut n = 0;

        for _ in 0..len {
            let bit = self.read_bit().map(u32::from)?;
            n = (n << 1) | bit;
        }

        Ok(n)
    }

    pub(crate) fn read_bit(&mut self) -> io::Result<u8> {
        if self.i >= 8 {
            let Some((b, rest)) = self.src.split_first() else {
                return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
            };

            self.buf = *b;
            self.src = rest;
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
        let src = [0b11001111, 0b01000000];
        let mut reader = BitReader::new(&src[..]);
        assert_eq!(reader.read_u32(4)?, 0b1100);
        assert_eq!(reader.read_u32(2)?, 0b11);
        assert_eq!(reader.read_u32(6)?, 0b110100);
        Ok(())
    }
}
