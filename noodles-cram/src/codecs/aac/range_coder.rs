#![allow(dead_code)]

use std::io::{self, Read};

use byteorder::ReadBytesExt;

#[derive(Debug)]
pub struct RangeCoder {
    range: u32,
    code: u32,
}

impl RangeCoder {
    pub fn range_decode_create<R>(&mut self, reader: &mut R) -> io::Result<()>
    where
        R: Read,
    {
        for _ in 0..=4 {
            let b = reader.read_u8().map(u32::from)?;
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
            let b = reader.read_u8().map(u32::from)?;
            self.code = (self.code << 8) | b;
        }

        Ok(())
    }
}

impl Default for RangeCoder {
    fn default() -> Self {
        Self {
            range: u32::MAX,
            code: 0,
        }
    }
}
