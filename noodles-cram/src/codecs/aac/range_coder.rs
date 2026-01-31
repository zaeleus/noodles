use std::io::{self, Write};

use crate::io::{reader::num::read_u8, writer::num::write_u8};

#[derive(Debug)]
pub struct RangeCoder {
    range: u32,
    code: u32,

    low: u32,
    carry: bool,
    cache: u32,
    ff_num: u32,
}

impl RangeCoder {
    pub fn new(src: &mut &[u8]) -> io::Result<Self> {
        // Discard first byte.
        read_u8(src)?;

        let code = read_u32_be(src)?;

        Ok(Self {
            code,
            ..Default::default()
        })
    }

    pub fn range_get_freq(&mut self, tot_freq: u32) -> u32 {
        self.range /= tot_freq;
        self.code / self.range
    }

    pub fn range_decode(&mut self, src: &mut &[u8], sym_low: u32, sym_freq: u32) -> io::Result<()> {
        self.code -= sym_low * self.range;
        self.range *= sym_freq;

        while self.range < (1 << 24) {
            self.range <<= 8;
            let b = read_u8(src).map(u32::from)?;
            self.code = (self.code << 8) | b;
        }

        Ok(())
    }

    pub fn range_encode<W>(
        &mut self,
        writer: &mut W,
        sym_low: u32,
        sym_freq: u32,
        tot_freq: u32,
    ) -> io::Result<()>
    where
        W: Write,
    {
        let old_low = self.low;

        self.range /= tot_freq;
        self.low = self.low.overflowing_add(sym_low * self.range).0;
        self.range *= sym_freq;

        if self.low < old_low {
            self.carry = true;
        }

        while self.range < (1 << 24) {
            self.range <<= 8;
            self.range_shift_low(writer)?;
        }

        Ok(())
    }

    fn range_shift_low<W>(&mut self, writer: &mut W) -> io::Result<()>
    where
        W: Write,
    {
        if self.low < 0xff000000 || self.carry {
            if !self.carry {
                let b = (self.cache & 0xff) as u8;
                write_u8(writer, b)?;

                while self.ff_num > 0 {
                    write_u8(writer, 0xff)?;
                    self.ff_num -= 1;
                }
            } else {
                let b = ((self.cache + 1) & 0xff) as u8;
                write_u8(writer, b)?;

                while self.ff_num > 0 {
                    write_u8(writer, 0x00)?;
                    self.ff_num -= 1;
                }
            }

            self.cache = self.low >> 24;
            self.carry = false;
        } else {
            self.ff_num += 1;
        }

        self.low <<= 8;

        Ok(())
    }

    pub fn range_encode_end<W>(&mut self, writer: &mut W) -> io::Result<()>
    where
        W: Write,
    {
        for _ in 0..=4 {
            self.range_shift_low(writer)?;
        }

        Ok(())
    }
}

impl Default for RangeCoder {
    fn default() -> Self {
        Self {
            range: u32::MAX,
            code: 0,

            low: 0,
            carry: false,
            cache: 0,
            ff_num: 0,
        }
    }
}

fn read_u32_be(src: &mut &[u8]) -> io::Result<u32> {
    let (buf, rest) = src
        .split_first_chunk()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Ok(u32::from_be_bytes(*buf))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() -> io::Result<()> {
        let src = [0x0d, 0x00, 0x00, 0x00, 0x08];
        let coder = RangeCoder::new(&mut &src[..])?;
        assert_eq!(coder.code, 8);

        assert!(matches!(
            RangeCoder::new(&mut &[][..]),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        Ok(())
    }
}
