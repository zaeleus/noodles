use std::{ops::Deref, slice};

const QUALITY_OFFSET: u8 = b'!';

#[derive(Debug)]
pub struct Quality<'a> {
    qual: &'a [u8],
}

impl<'a> Quality<'a> {
    pub fn new(qual: &[u8]) -> Quality {
        Quality { qual }
    }

    pub fn chars(&self) -> Chars<slice::Iter<u8>> {
        Chars {
            chars: self.qual.iter(),
        }
    }
}

impl<'a> Deref for Quality<'a> {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.qual
    }
}

pub struct Chars<I> {
    chars: I,
}

impl<'a, I: Iterator<Item = &'a u8>> Iterator for Chars<I> {
    type Item = char;

    fn next(&mut self) -> Option<char> {
        self.chars.next().map(|b| (b + QUALITY_OFFSET) as char)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.chars.size_hint()
    }
}

impl<'a, I: Iterator<Item = &'a u8> + DoubleEndedIterator> DoubleEndedIterator for Chars<I> {
    fn next_back(&mut self) -> Option<char> {
        self.chars.next_back().map(|b| (b + QUALITY_OFFSET) as char)
    }
}

#[cfg(test)]
mod tests {
    use super::{Quality, QUALITY_OFFSET};

    #[test]
    fn test_chars() {
        let data: Vec<_> = b"><>=@>;".iter().map(|b| b - QUALITY_OFFSET).collect();
        let quality = Quality::new(&data);
        let actual: Vec<char> = quality.chars().collect();
        assert_eq!(actual, vec!['>', '<', '>', '=', '@', '>', ';']);
    }
}
