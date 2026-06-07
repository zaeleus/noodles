use std::{io, num::NonZero};

use bstr::{BStr, BString};
use noodles_core::region::Interval;

/// A FASTA index record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Record {
    name: BString,
    length: u64,
    position: u64,
    line_base_count: NonZero<u64>,
    line_width: NonZero<u64>,
}

impl Record {
    /// Creates a FASTA index record.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZero;
    ///
    /// use noodles_fasta::fai;
    ///
    /// let line_base_count = const { NonZero::new(80).unwrap() };
    /// let line_width = const { NonZero::new(81).unwrap() };
    /// let record = fai::Record::new("sq0", 8, 4, line_base_count, line_width);
    /// ```
    pub fn new<N>(
        name: N,
        length: u64,
        position: u64,
        line_base_count: NonZero<u64>,
        line_width: NonZero<u64>,
    ) -> Self
    where
        N: Into<BString>,
    {
        Self {
            name: name.into(),
            length,
            position,
            line_base_count,
            line_width,
        }
    }

    /// Returns the record name.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZero;
    ///
    /// use noodles_fasta::fai;
    ///
    /// let line_base_count = const { NonZero::new(80).unwrap() };
    /// let line_width = const { NonZero::new(81).unwrap() };
    /// let record = fai::Record::new("sq0", 8, 4, line_base_count, line_width);
    /// assert_eq!(record.name(), b"sq0");
    /// ```
    pub fn name(&self) -> &BStr {
        self.name.as_ref()
    }

    /// Returns the length of the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZero;
    ///
    /// use noodles_fasta::fai;
    ///
    /// let line_base_count = const { NonZero::new(80).unwrap() };
    /// let line_width = const { NonZero::new(81).unwrap() };
    /// let record = fai::Record::new("sq0", 8, 4, line_base_count, line_width);
    /// assert_eq!(record.length(), 8);
    /// ```
    pub fn length(&self) -> u64 {
        self.length
    }

    /// Returns the position of the start of the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZero;
    ///
    /// use noodles_fasta::fai;
    ///
    /// let line_base_count = const { NonZero::new(80).unwrap() };
    /// let line_width = const { NonZero::new(81).unwrap() };
    /// let record = fai::Record::new("sq0", 8, 4, line_base_count, line_width);
    /// assert_eq!(record.position(), 4);
    /// ```
    pub fn position(&self) -> u64 {
        self.position
    }

    /// Returns the offset from the start.
    #[deprecated(since = "0.63.0", note = "Use `Record::position` instead.")]
    pub fn offset(&self) -> u64 {
        self.position()
    }

    /// Returns the number of bases in a line.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZero;
    ///
    /// use noodles_fasta::fai;
    ///
    /// let line_base_count = const { NonZero::new(80).unwrap() };
    /// let line_width = const { NonZero::new(81).unwrap() };
    /// let record = fai::Record::new("sq0", 8, 4, line_base_count, line_width);
    /// assert_eq!(record.line_base_count(), line_base_count);
    /// ```
    pub fn line_base_count(&self) -> NonZero<u64> {
        self.line_base_count
    }

    /// Returns the number of bases in a line.
    #[deprecated(since = "0.63.0", note = "Use `Record::line_base_count` instead.")]
    pub fn line_bases(&self) -> u64 {
        self.line_base_count.get()
    }

    /// Returns the number of characters in a line.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZero;
    ///
    /// use noodles_fasta::fai;
    ///
    /// let line_base_count = const { NonZero::new(80).unwrap() };
    /// let line_width = const { NonZero::new(81).unwrap() };
    /// let record = fai::Record::new("sq0", 8, 4, line_base_count, line_width);
    /// assert_eq!(record.line_width(), line_width);
    /// ```
    pub fn line_width(&self) -> NonZero<u64> {
        self.line_width
    }

    /// Returns the start position of the given interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZero;
    ///
    /// use noodles_core::{region::Interval, Position};
    /// use noodles_fasta::fai;
    ///
    /// let line_base_count = const { NonZero::new(80).unwrap() };
    /// let line_width = const { NonZero::new(81).unwrap() };
    /// let record = fai::Record::new("sq0", 8, 4, line_base_count, line_width);
    /// let interval = Interval::from(..);
    ///
    /// assert_eq!(record.query(interval)?, 4);
    /// Ok::<_, std::io::Error>(())
    /// ```
    pub fn query(&self, interval: Interval) -> io::Result<u64> {
        let start = interval
            .start()
            .map(|position| usize::from(position) - 1)
            .unwrap_or_default();

        let start =
            u64::try_from(start).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let line_base_count = self.line_base_count.get();
        let line_width = self.line_width.get();
        let pos = self.position() + start / line_base_count * line_width + start % line_base_count;

        Ok(pos)
    }
}

impl Default for Record {
    fn default() -> Self {
        Self {
            name: BString::default(),
            length: 0,
            position: 0,
            line_base_count: NonZero::<u64>::MIN,
            line_width: NonZero::<u64>::MIN,
        }
    }
}
