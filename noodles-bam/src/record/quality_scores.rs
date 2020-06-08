use std::{ops::Deref, slice};

const QUALITY_OFFSET: u8 = b'!';

#[derive(Debug)]
pub struct QualityScores<'a> {
    qual: &'a [u8],
}

impl<'a> QualityScores<'a> {
    pub fn new(qual: &'a [u8]) -> Self {
        Self { qual }
    }

    pub fn chars(&self) -> Chars<slice::Iter<'_, u8>> {
        Chars {
            chars: self.qual.iter(),
        }
    }
}

impl<'a> Deref for QualityScores<'a> {
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
        self.chars.next().map(|&b| byte_to_char(b))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.chars.size_hint()
    }
}

impl<'a, I: Iterator<Item = &'a u8> + DoubleEndedIterator> DoubleEndedIterator for Chars<I> {
    fn next_back(&mut self) -> Option<char> {
        self.chars.next_back().map(|&b| byte_to_char(b))
    }
}

fn byte_to_char(b: u8) -> char {
    (b + QUALITY_OFFSET) as char
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chars() {
        let data: Vec<_> = b"><>=@>;".iter().map(|b| b - QUALITY_OFFSET).collect();
        let quality = QualityScores::new(&data);
        let actual: Vec<char> = quality.chars().collect();
        assert_eq!(actual, vec!['>', '<', '>', '=', '@', '>', ';']);
    }
}
