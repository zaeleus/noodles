use std::io::{self, Write};

use super::RangeCoder;

#[derive(Clone, Debug)]
pub struct Model {
    total_freq: u32,
    symbols: Vec<u8>,
    frequencies: Vec<u32>,
}

impl Model {
    pub fn new(max_sym: u8) -> Self {
        let num_sym = usize::from(max_sym) + 1;

        let mut symbols = Vec::with_capacity(num_sym);

        for i in 0..=max_sym {
            symbols.push(i);
        }

        let frequencies = vec![1; num_sym];

        Self {
            total_freq: u32::from(max_sym) + 1,
            symbols,
            frequencies,
        }
    }

    pub fn decode(&mut self, src: &mut &[u8], range_coder: &mut RangeCoder) -> io::Result<u8> {
        let freq = range_coder.range_get_freq(self.total_freq);

        let mut acc = 0;
        let mut x = 0;

        while acc + self.frequencies[x] <= freq {
            acc += self.frequencies[x];
            x += 1;
        }

        range_coder.range_decode(src, acc, self.frequencies[x])?;

        self.frequencies[x] += 16;
        self.total_freq += 16;

        if self.total_freq > (1 << 16) - 17 {
            self.renormalize();
        }

        let sym = self.symbols[x];

        if x > 0 && self.frequencies[x] > self.frequencies[x - 1] {
            self.frequencies.swap(x, x - 1);
            self.symbols.swap(x, x - 1);
        }

        Ok(sym)
    }

    pub fn encode<W>(
        &mut self,
        writer: &mut W,
        range_coder: &mut RangeCoder,
        sym: u8,
    ) -> io::Result<()>
    where
        W: Write,
    {
        let mut acc = 0;
        let mut x = 0;

        while self.symbols[x] != sym {
            acc += self.frequencies[x];
            x += 1;
        }

        let sym_freq = self.frequencies[x];
        range_coder.range_encode(writer, acc, sym_freq, self.total_freq)?;

        self.frequencies[x] += 16;
        self.total_freq += 16;

        if self.total_freq > (1 << 16) - 17 {
            self.renormalize();
        }

        if x > 0 && self.frequencies[x] > self.frequencies[x - 1] {
            self.frequencies.swap(x, x - 1);
            self.symbols.swap(x, x - 1);
        }

        Ok(())
    }

    fn renormalize(&mut self) {
        let mut total_freq = 0;

        for freq in &mut self.frequencies {
            *freq -= *freq / 2;
            total_freq += *freq;
        }

        self.total_freq = total_freq;
    }
}
