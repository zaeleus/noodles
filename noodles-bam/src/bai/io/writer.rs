mod index;

use std::io::{self, Write};

use self::index::write_index;
use crate::bai::Index;

/// A BAM index (BAI) writer.
///
/// # Examples
///
/// ```no_run
/// # use std::{fs::File, io};
/// use noodles_bam::bai;
/// use noodles_csi as csi;
///
/// let index = bai::Index::default();
///
/// let mut writer = File::create("sample.bam.bai").map(bai::io::Writer::new)?;
/// writer.write_index(&index)?;
/// # Ok::<(), io::Error>(())
/// ```
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::bai;
    /// let writer = bai::io::Writer::new(io::sink());
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
    /// use noodles_bam::bai;
    /// let mut writer = bai::io::Writer::new(io::sink());
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
    /// use noodles_bam::bai;
    /// let writer = bai::io::Writer::new(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a BAI index writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// let writer = bai::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Writes a BAM index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::bai;
    ///
    /// let index = bai::Index::default();
    ///
    /// let mut writer = bai::io::Writer::new(Vec::new());
    /// writer.write_index(&index)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_index(&mut self, index: &Index) -> io::Result<()> {
        write_index(&mut self.inner, index)
    }
}
