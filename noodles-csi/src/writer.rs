pub(crate) mod index;

use std::io::{self, Write};

use noodles_bgzf as bgzf;

use self::index::write_index;
use super::{index::reference_sequence::index::BinnedIndex, Index};

/// A CSI writer.
pub struct Writer<W>
where
    W: Write,
{
    inner: bgzf::Writer<W>,
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
    /// let writer = csi::Writer::new(Vec::new());
    /// ```
    pub fn new(writer: W) -> Self {
        Self {
            inner: bgzf::Writer::new(writer),
        }
    }

    /// Writes a coordinate-sorted index (CSI).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// let mut writer = csi::Writer::new(Vec::new());
    /// writer.write_index(&index)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_index(&mut self, index: &Index<BinnedIndex>) -> io::Result<()> {
        write_index(&mut self.inner, index)
    }
}
