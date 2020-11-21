use std::{ops::Deref, slice};

const QUALITY_OFFSET: u8 = b'!';

/// BAM record quality scores.
#[derive(Debug)]
pub struct QualityScores<'a> {
    qual: &'a [u8],
}

impl<'a> QualityScores<'a> {
    /// Creates quality scores from raw quality scores data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    ///
    /// // NDLS
    /// let data = [45, 35, 43, 50];
    /// let quality_scores = QualityScores::new(&data);
    /// ```
    pub fn new(qual: &'a [u8]) -> Self {
        Self { qual }
    }

    /// Returns an iterator over quality scores as offset printable ASCII characters.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::QualityScores;
    ///
    /// // NDLS
    /// let data = [45, 35, 43, 50];
    /// let quality_scores = QualityScores::new(&data);
    ///
    /// let mut chars = quality_scores.chars();
    ///
    /// assert_eq!(chars.next(), Some('N'));
    /// assert_eq!(chars.next(), Some('D'));
    /// assert_eq!(chars.next(), Some('L'));
    /// assert_eq!(chars.next(), Some('S'));
    /// assert_eq!(chars.next(), None);
    /// ```
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

/// An iterator over quality scores as offset printable ASCII characters.
///
/// This is created by calling [`QualityScores::chars`].
pub struct Chars<I> {
    chars: I,
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
    use super::*;

    #[test]
    fn test_chars() {
        let data: Vec<_> = b"><>=@>;".iter().map(|b| b - QUALITY_OFFSET).collect();
        let quality = QualityScores::new(&data);
        let actual: Vec<char> = quality.chars().collect();
        assert_eq!(actual, vec!['>', '<', '>', '=', '@', '>', ';']);
    }
}
