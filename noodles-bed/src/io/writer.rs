mod record;

use std::io::{self, Write};

use self::record::{write_record_3, write_record_4, write_record_5, write_record_6};
use crate::Record;

/// A BED writer.
pub struct Writer<W, const N: usize> {
    inner: W,
}

impl<W, const N: usize> Writer<W, N> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let writer = bed::io::Writer::<_, 3>::new(Vec::new());
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
    /// let mut writer = bed::io::Writer::<_, 3>::new(Vec::new());
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
    /// let writer = bed::io::Writer::<_, 3>::new(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }
}

impl<W, const N: usize> Writer<W, N>
where
    W: Write,
{
    /// Creates a BED writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let writer = bed::io::Writer::<_, 3>::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }
}

impl<W> Writer<W, 3>
where
    W: Write,
{
    /// Writes a record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<_, 3>::new(Vec::new());
    /// let record = bed::Record::default();
    /// writer.write_record(&record)?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record<3>) -> io::Result<()> {
        self.write_feature_record(record)
    }

    /// Writes a feature record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<_, 3>::new(Vec::new());
    /// let record = bed::Record::default();
    /// writer.write_feature_record(&record)?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn write_feature_record<R>(&mut self, record: &R) -> io::Result<()>
    where
        R: crate::feature::Record<3>,
    {
        write_record_3(&mut self.inner, record)
    }
}

impl<W> Writer<W, 4>
where
    W: Write,
{
    /// Writes a record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<_, 4>::new(Vec::new());
    /// let record = bed::Record::default();
    /// writer.write_record(&record)?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record<4>) -> io::Result<()> {
        self.write_feature_record(record)
    }

    /// Writes a feature record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<_, 4>::new(Vec::new());
    /// let record = bed::Record::default();
    /// writer.write_feature_record(&record)?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn write_feature_record<R>(&mut self, record: &R) -> io::Result<()>
    where
        R: crate::feature::Record<4>,
    {
        write_record_4(&mut self.inner, record)
    }
}

impl<W> Writer<W, 5>
where
    W: Write,
{
    /// Writes a record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<_, 5>::new(Vec::new());
    /// let record = bed::Record::default();
    /// writer.write_record(&record)?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record<5>) -> io::Result<()> {
        self.write_feature_record(record)
    }

    /// Writes a feature record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<_, 5>::new(Vec::new());
    /// let record = bed::Record::default();
    /// writer.write_feature_record(&record)?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn write_feature_record<R>(&mut self, record: &R) -> io::Result<()>
    where
        R: crate::feature::Record<5>,
    {
        write_record_5(&mut self.inner, record)
    }
}

impl<W> Writer<W, 6>
where
    W: Write,
{
    /// Writes a record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<_, 6>::new(Vec::new());
    /// let record = bed::Record::default();
    /// writer.write_record(&record)?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record<6>) -> io::Result<()> {
        self.write_feature_record(record)
    }

    /// Writes a feature record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<_, 6>::new(Vec::new());
    /// let record = bed::Record::default();
    /// writer.write_feature_record(&record)?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn write_feature_record<R>(&mut self, record: &R) -> io::Result<()>
    where
        R: crate::feature::Record<6>,
    {
        write_record_6(&mut self.inner, record)
    }
}
