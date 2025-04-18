use std::io;

/// SAM record quality scores.
#[derive(Debug, Eq, PartialEq)]
pub struct QualityScores<'a>(&'a [u8]);

impl<'a> QualityScores<'a> {
    /// Creates SAM record quality scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::QualityScores;
    /// let quality_scores = QualityScores::new(b"NDLS");
    /// ```
    pub fn new(buf: &'a [u8]) -> Self {
        Self(buf)
    }

    /// Returns whether there are any scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::QualityScores;
    ///
    /// let quality_scores = QualityScores::new(b"");
    /// assert!(quality_scores.is_empty());
    ///
    /// let quality_scores = QualityScores::new(b"NDLS");
    /// assert!(!quality_scores.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of scores.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::record::QualityScores;
    /// let quality_scores = QualityScores::new(b"NDLS");
    /// assert_eq!(quality_scores.len(), 4);
    /// ```
    pub fn len(&self) -> usize {
        self.0.len()
    }
}

impl crate::alignment::record::QualityScores for QualityScores<'_> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u8>> + '_> {
        const OFFSET: u8 = b'!';

        Box::new(self.as_ref().iter().map(|&b| {
            b.checked_sub(OFFSET)
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid score"))
        }))
    }
}

impl AsRef<[u8]> for QualityScores<'_> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use crate::alignment::record::QualityScores as _;

    use super::*;

    #[test]
    fn test_iter() -> io::Result<()> {
        let quality_scores = QualityScores::new(b"NDLS");

        assert_eq!(
            quality_scores.iter().collect::<io::Result<Vec<_>>>()?,
            [45, 35, 43, 50]
        );

        let quality_scores = QualityScores::new(&[0x00]);

        assert!(matches!(
            quality_scores.iter().collect::<io::Result<Vec<_>>>(),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
