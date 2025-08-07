pub(crate) mod index;
mod num;

use std::io::{self, Write};

use noodles_bgzf as bgzf;

use self::index::write_index;
use crate::Index;

/// A CSI writer.
pub struct Writer<W>
where
    W: Write,
{
    inner: bgzf::io::Writer<W>,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a CSI writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let writer = csi::io::Writer::new(Vec::new());
    /// ```
    pub fn new(writer: W) -> Self {
        Self {
            inner: bgzf::io::Writer::new(writer),
        }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    /// let writer = csi::io::Writer::new(io::sink());
    /// let _inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &bgzf::io::Writer<W> {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    /// let mut writer = csi::io::Writer::new(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut bgzf::io::Writer<W> {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    /// let writer = csi::io::Writer::new(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> bgzf::io::Writer<W> {
        self.inner
    }

    /// Writes a coordinate-sorted index (CSI).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// let mut writer = csi::io::Writer::new(Vec::new());
    /// writer.write_index(&index)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_index(&mut self, index: &Index) -> io::Result<()> {
        write_index(&mut self.inner, index)
    }
}
