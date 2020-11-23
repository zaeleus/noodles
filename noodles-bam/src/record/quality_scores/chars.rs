const QUALITY_OFFSET: u8 = b'!';

/// An iterator over quality scores as offset printable ASCII characters.
///
/// This is created by calling [`super::QualityScores::chars`].
pub struct Chars<I> {
    chars: I,
}

impl<I> Chars<I> {
    pub(crate) fn new(iter: I) -> Self {
        Self { chars: iter }
    }
}

impl<'a, I> Iterator for Chars<I>
where
    I: Iterator<Item = &'a u8>,
{
    type Item = char;

    fn next(&mut self) -> Option<char> {
        self.chars.next().map(|&b| byte_to_char(b))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.chars.size_hint()
    }
}

impl<'a, I> DoubleEndedIterator for Chars<I>
where
    I: Iterator<Item = &'a u8> + DoubleEndedIterator,
{
    fn next_back(&mut self) -> Option<char> {
        self.chars.next_back().map(|&b| byte_to_char(b))
    }
}

fn byte_to_char(b: u8) -> char {
    (b + QUALITY_OFFSET) as char
}

#[cfg(test)]
mod tests {
    use crate::record::QualityScores;

    use super::*;

    #[test]
    fn test_chars() {
        let data: Vec<_> = b"><>=@>;".iter().map(|b| b - QUALITY_OFFSET).collect();
        let quality = QualityScores::new(&data);
        let actual: Vec<char> = quality.chars().collect();
        assert_eq!(actual, vec!['>', '<', '>', '=', '@', '>', ';']);
    }
}
