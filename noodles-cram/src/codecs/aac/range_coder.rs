use std::io::{self, Read, Write};

use byteorder::WriteBytesExt;

use crate::io::reader::num::read_u8;

#[derive(Debug)]
pub struct RangeCoder {
    range: u32,
    code: u32,

    low: u32,
    carry: u32,
    cache: u32,
    ff_num: u32,
}

impl RangeCoder {
    pub fn range_decode_create<R>(&mut self, reader: &mut R) -> io::Result<()>
    where
        R: Read,
    {
        for _ in 0..=4 {
            let b = read_u8(reader).map(u32::from)?;
            self.code = (self.code << 8) | b;
        }

        self.code &= u32::MAX;

        Ok(())
    }

    pub fn range_get_freq(&mut self, tot_freq: u32) -> u32 {
        self.range /= tot_freq;
        self.code / self.range
    }

    pub fn range_decode<R>(&mut self, reader: &mut R, sym_low: u32, sym_freq: u32) -> io::Result<()>
    where
        R: Read,
    {
        self.code -= sym_low * self.range;
        self.range *= sym_freq;

        while self.range < (1 << 24) {
            self.range <<= 8;
            let b = read_u8(reader).map(u32::from)?;
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
            self.carry = 1;
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
        if self.low < 0xff000000 || self.carry != 0 {
            if self.carry == 0 {
                let b = (self.cache & 0xff) as u8;
                writer.write_u8(b)?;

                while self.ff_num > 0 {
                    writer.write_u8(0xff)?;
                    self.ff_num -= 1;
                }
            } else {
                let b = ((self.cache + 1) & 0xff) as u8;
                writer.write_u8(b)?;

                while self.ff_num > 0 {
                    writer.write_u8(0x00)?;
                    self.ff_num -= 1;
                }
            }

            self.cache = self.low >> 24;
            self.carry = 0;
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
            carry: 0,
            cache: 0,
            ff_num: 0,
        }
    }
}
