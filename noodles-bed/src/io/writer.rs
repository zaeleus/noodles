//! BED writer.

mod builder;
mod record;

use std::io::{self, Write};

pub use self::builder::Builder;
use self::record::{write_record_3, write_record_4, write_record_5, write_record_6};
use crate::Record;

/// A BED writer.
pub struct Writer<const N: usize, W> {
    inner: W,
}

impl<const N: usize, W> Writer<N, W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed as bed;
    /// let writer = bed::io::Writer::<3, _>::new(io::sink());
    /// let _inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<3, _>::new(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed as bed;
    /// let writer = bed::io::Writer::<3, _>::new(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }
}

impl<const N: usize, W> Writer<N, W>
where
    W: Write,
{
    /// Creates a BED writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed as bed;
    /// let writer = bed::io::Writer::<3, _>::new(io::sink());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }
}

impl<W> Writer<3, W>
where
    W: Write,
{
    /// Writes a record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<3, _>::new(io::sink());
    /// let record = bed::Record::default();
    /// writer.write_record(&record)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record<3>) -> io::Result<()> {
        self.write_feature_record(record)
    }

    /// Writes a feature record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<3, _>::new(io::sink());
    /// let record = bed::Record::default();
    /// writer.write_feature_record(&record)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_feature_record<R>(&mut self, record: &R) -> io::Result<()>
    where
        R: crate::feature::Record<3>,
    {
        write_record_3(&mut self.inner, record)
    }
}

impl<W> Writer<4, W>
where
    W: Write,
{
    /// Writes a record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<4, _>::new(io::sink());
    /// let record = bed::Record::default();
    /// writer.write_record(&record)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record<4>) -> io::Result<()> {
        self.write_feature_record(record)
    }

    /// Writes a feature record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<4, _>::new(io::sink());
    /// let record = bed::Record::default();
    /// writer.write_feature_record(&record)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_feature_record<R>(&mut self, record: &R) -> io::Result<()>
    where
        R: crate::feature::Record<4>,
    {
        write_record_4(&mut self.inner, record)
    }
}

impl<W> Writer<5, W>
where
    W: Write,
{
    /// Writes a record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<5, _>::new(io::sink());
    /// let record = bed::Record::default();
    /// writer.write_record(&record)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record<5>) -> io::Result<()> {
        self.write_feature_record(record)
    }

    /// Writes a feature record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<5, _>::new(io::sink());
    /// let record = bed::Record::default();
    /// writer.write_feature_record(&record)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_feature_record<R>(&mut self, record: &R) -> io::Result<()>
    where
        R: crate::feature::Record<5>,
    {
        write_record_5(&mut self.inner, record)
    }
}

impl<W> Writer<6, W>
where
    W: Write,
{
    /// Writes a record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<6, _>::new(io::sink());
    /// let record = bed::Record::default();
    /// writer.write_record(&record)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record<6>) -> io::Result<()> {
        self.write_feature_record(record)
    }

    /// Writes a feature record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bed as bed;
    /// let mut writer = bed::io::Writer::<6, _>::new(io::sink());
    /// let record = bed::Record::default();
    /// writer.write_feature_record(&record)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_feature_record<R>(&mut self, record: &R) -> io::Result<()>
    where
        R: crate::feature::Record<6>,
    {
        write_record_6(&mut self.inner, record)
    }
}
