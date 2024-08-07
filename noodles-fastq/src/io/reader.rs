//! FASTQ reader.

pub(crate) mod record;
mod records;

pub use self::records::Records;

use std::io::{self, BufRead};

use self::record::read_record;
use crate::Record;

/// A FASTQ reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fastq as fastq;
    /// let reader = fastq::io::Reader::new(io::empty());
    /// let _inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fastq as fastq;
    /// let mut reader = fastq::io::Reader::new(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Unwraps and returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fastq as fastq;
    /// let reader = fastq::io::Reader::new(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> Reader<R>
where
    R: BufRead,
{
    /// Creates a FASTQ reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// let data = b"@r0\nATCG\n+\nNDLS\n";
    /// let reader = fastq::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a FASTQ record.
    ///
    /// This reads from the underlying stream until four lines are read: the read name, the
    /// sequence, the plus line, and the quality scores. Each line omits the trailing newline.
    ///
    /// The stream is expected to be at the start of a record.
    ///
    /// If successful, the number of bytes read is returned. If the number of bytes read is 0, the
    /// stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fastq as fastq;
    ///
    /// let data = b"@r0\nATCG\n+\nNDLS\n";
    /// let mut reader = fastq::io::Reader::new(&data[..]);
    ///
    /// let mut record = fastq::Record::default();
    /// reader.read_record(&mut record)?;
    ///
    /// assert_eq!(record.name(), &b"r0"[..]);
    /// assert_eq!(record.sequence(), b"ATCG");
    /// assert_eq!(record.quality_scores(), b"NDLS");
    /// Ok::<(), io::Error>(())
    /// ```
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        read_record(&mut self.inner, record)
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// The stream is expected to be at the start of a record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fastq::{self as fastq, record::Definition};
    ///
    /// let data = b"@r0\nATCG\n+\nNDLS\n";
    /// let mut reader = fastq::io::Reader::new(&data[..]);
    ///
    /// let mut records = reader.records();
    ///
    /// assert_eq!(
    ///     records.next().transpose()?,
    ///     Some(fastq::Record::new(Definition::new("r0", ""), "ATCG", "NDLS")
    /// ));
    ///
    /// assert!(records.next().is_none());
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn records(&mut self) -> Records<'_, R> {
        Records::new(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::Definition;

    #[test]
    fn test_read_record() -> io::Result<()> {
        let data = b"\
@noodles:1/1
AGCT
+
abcd
@noodles:2/1
TCGA
+noodles:2/1
dcba
";

        let mut reader = &data[..];
        let mut record = Record::default();

        read_record(&mut reader, &mut record)?;
        let expected = Record::new(Definition::new("noodles:1/1", ""), "AGCT", "abcd");
        assert_eq!(record, expected);

        read_record(&mut reader, &mut record)?;
        let expected = Record::new(Definition::new("noodles:2/1", ""), "TCGA", "dcba");
        assert_eq!(record, expected);

        let n = read_record(&mut reader, &mut record)?;
        assert_eq!(n, 0);

        Ok(())
    }
}
