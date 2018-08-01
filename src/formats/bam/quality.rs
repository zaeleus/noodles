use std::slice;

const QUALITY_OFFSET: u8 = b'!';

#[derive(Debug)]
pub struct Quality {
    qual: Vec<u8>,
}

impl Quality {
    pub fn new(qual: Vec<u8>) -> Quality {
        Quality { qual }
    }

    pub fn len(&self) -> usize {
        self.qual.len()
    }

    pub fn chars(&self) -> Chars<slice::Iter<u8>> {
        Chars { chars: self.qual.iter() }
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.qual
    }
}

pub struct Chars<I> {
    chars: I,
}

impl<'a, I: Iterator<Item=&'a u8>> Iterator for Chars<I> {
    type Item = char;

    fn next(&mut self) -> Option<char> {
        self.chars.next().map(|b| (b + QUALITY_OFFSET) as char)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.chars.size_hint()
    }
}

impl<'a, I: Iterator<Item=&'a u8> + DoubleEndedIterator> DoubleEndedIterator for Chars<I> {
    fn next_back(&mut self) -> Option<char> {
        self.chars.next_back().map(|b| (b + QUALITY_OFFSET) as char)
    }
}
