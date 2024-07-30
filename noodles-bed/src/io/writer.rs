mod record;

use std::io::{self, Write};

use self::record::{write_record_3, write_record_4, write_record_5, write_record_6};
use crate::Record;

/// A BED writer.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let writer = bed::io::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::new(Vec::new());
    /// assert!(writer.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let writer = bed::io::Writer::new(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a BED writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let writer = bed::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Writes a record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::new(Vec::new());
    /// let record = bed::Record::<3>::default();
    /// writer.write_record(&record)?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn write_record<const N: usize>(&mut self, record: &Record<N>) -> io::Result<()>
    where
        Record<N>: crate::feature::Record<N>,
    {
        self.write_feature_record(record)
    }

    /// Writes a feature record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::new(Vec::new());
    /// let record = bed::Record::<3>::default();
    /// writer.write_feature_record(&record)?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn write_feature_record<R, const N: usize>(&mut self, record: &R) -> io::Result<()>
    where
        R: crate::feature::Record<N>,
    {
        match record.standard_field_count() {
            3 => write_record_3(&mut self.inner, record),
            4 => write_record_4(&mut self.inner, record),
            5 => write_record_5(&mut self.inner, record),
            6 => write_record_6(&mut self.inner, record),
            _ => todo!(),
        }
    }
}
