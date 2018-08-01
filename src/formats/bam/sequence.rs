use std::ops::Index;

static SEQ_CHARS: &[char] = &[
    '=', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
    'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N',
];

#[derive(Debug)]
pub struct Sequence {
    seq: Vec<u8>,
    len: usize,
}

impl Sequence {
    pub fn new(seq: Vec<u8>, len: usize) -> Sequence {
        Sequence { seq, len }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn symbols(&self) -> Symbols {
        Symbols {
            sequence: self,
            head: 0,
            tail: self.len() - 1,
            remaining: self.len(),
        }
    }
}

pub struct Symbols<'a> {
    sequence: &'a Sequence,
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

impl<I: Iterator<Item=char>> Complement<I> {
    pub fn new(iter: I) -> Complement<I> {
        Complement { iter }
    }
}

impl<I: Iterator<Item=char>> Iterator for Complement<I> {
    type Item = char;

    fn next(&mut self) -> Option<char> {
        self.iter.next().map(|b| complement(b))
    }
}

impl Index<usize> for Sequence {
    type Output = char;

    fn index(&self, i: usize) -> &char {
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
