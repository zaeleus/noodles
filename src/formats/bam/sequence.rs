use std::ops::{Deref, Index};

static SEQ_CHARS: &[char] = &[
    '=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N',
];

#[derive(Debug)]
pub struct Sequence<'a> {
    seq: &'a [u8],
    n_chars: usize,
}

impl<'a> Sequence<'a> {
    pub fn new(seq: &[u8], n_chars: usize) -> Sequence {
        Sequence { seq, n_chars }
    }

    pub fn n_chars(&self) -> usize {
        self.n_chars
    }

    pub fn symbols(&self) -> Symbols {
        Symbols {
            sequence: self,
            head: 0,
            tail: self.n_chars - 1,
            remaining: self.n_chars,
        }
    }
}

impl<'a> Deref for Sequence<'a> {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.seq
    }
}

pub struct Symbols<'a> {
    sequence: &'a Sequence<'a>,
    head: usize,
    tail: usize,
    remaining: usize,
}

impl<'a> Iterator for Symbols<'a> {
    type Item = char;

    fn next(&mut self) -> Option<char> {
        if self.remaining == 0 {
            return None;
        }

        let symbol = self.sequence[self.head];
        self.head += 1;
        self.remaining -= 1;
        Some(symbol)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl<'a> DoubleEndedIterator for Symbols<'a> {
    fn next_back(&mut self) -> Option<char> {
        if self.remaining == 0 {
            return None;
        }

        let symbol = self.sequence[self.tail];
        self.tail -= 1;
        self.remaining -= 1;
        Some(symbol)
    }
}

pub struct Complement<I> {
    iter: I,
}

impl<I: Iterator<Item = char>> Complement<I> {
    pub fn new(iter: I) -> Complement<I> {
        Complement { iter }
    }
}

impl<I: Iterator<Item = char>> Iterator for Complement<I> {
    type Item = char;

    fn next(&mut self) -> Option<char> {
        self.iter.next().map(|b| complement(b))
    }
}

impl<'a> Index<usize> for Sequence<'a> {
    type Output = char;

    fn index(&self, i: usize) -> &Self::Output {
        let j = i / 2;
        let b = self.seq[j];

        let k = if i % 2 == 0 {
            (b & 0xf0) >> 4
        } else {
            b & 0x0f
        };

        &SEQ_CHARS[k as usize]
    }
}

// https://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation
fn complement(s: char) -> char {
    match s {
        '=' => '=',
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        'W' => 'W',
        'S' => 'S',
        'M' => 'K',
        'K' => 'M',
        'R' => 'Y',
        'Y' => 'R',
        'B' => 'V',
        'D' => 'H',
        'H' => 'D',
        'V' => 'B',
        'N' => 'N',
        _ => panic!("invalid symbol '{}'", s),
    }
}

#[cfg(test)]
mod tests {
    use super::complement;

    #[test]
    fn test_complement() {
        assert_eq!(complement('='), '=');
        assert_eq!(complement('A'), 'T');
        assert_eq!(complement('C'), 'G');
        assert_eq!(complement('G'), 'C');
        assert_eq!(complement('T'), 'A');
        assert_eq!(complement('W'), 'W');
        assert_eq!(complement('S'), 'S');
        assert_eq!(complement('M'), 'K');
        assert_eq!(complement('K'), 'M');
        assert_eq!(complement('R'), 'Y');
        assert_eq!(complement('Y'), 'R');
        assert_eq!(complement('B'), 'V');
        assert_eq!(complement('D'), 'H');
        assert_eq!(complement('H'), 'D');
        assert_eq!(complement('V'), 'B');
        assert_eq!(complement('N'), 'N');
    }

    #[test]
    #[should_panic]
    fn test_complement_with_invalid_symbol() {
        complement('?');
    }
}
