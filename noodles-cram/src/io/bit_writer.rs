use std::io;

#[derive(Debug, Default)]
pub struct BitWriter {
    dst: Vec<u8>,
    buf: u8,
    i: usize,
}

impl BitWriter {
    pub fn flush(&mut self) -> io::Result<()> {
        if self.i > 0 {
            self.write_u32(0, 8 - self.i)
        } else {
            Ok(())
        }
    }

    pub fn finish(mut self) -> io::Result<Vec<u8>> {
        self.flush()?;
        Ok(self.dst)
    }

    pub fn write_u32(&mut self, value: u32, len: usize) -> io::Result<()> {
        if len == 0 {
            return Ok(());
        } else if len > 32 {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }

        let mut mask = 0x01 << (len - 1);

        for _ in 0..len {
            self.write_bit(value & mask != 0);
            mask >>= 1;
        }

        Ok(())
    }

    fn write_bit(&mut self, is_set: bool) {
        if is_set {
            self.buf |= 0x01 << (8 - self.i - 1);
        }

        self.i += 1;

        if self.i == 8 {
            self.dst.push(self.buf);
            self.buf = 0;
            self.i = 0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_u32() -> io::Result<()> {
        let mut writer = BitWriter::default();

        writer.write_u32(0x0c, 4)?;
        writer.write_u32(0x03, 2)?;
        writer.write_u32(0x34, 6)?;
        writer.flush()?;

        let expected = [0b11001111, 0b01000000];
        assert_eq!(writer.dst, expected);

        Ok(())
    }

    #[test]
    fn test_write_u32_with_0_len() -> io::Result<()> {
        let mut writer = BitWriter::default();
        writer.write_u32(0xff, 0)?;
        assert!(writer.dst.is_empty());
        Ok(())
    }

    #[test]
    fn test_write_u32_with_length_greater_than_32_bits() {
        let mut writer = BitWriter::default();
        assert!(writer.write_u32(0xff, 33).is_err());
    }
}
