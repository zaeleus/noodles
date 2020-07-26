use std::io::{self, Write};

use super::MAGIC_NUMBER;

// [major, minor]
const FILE_DEFINITION_FORMAT: [u8; 2] = [3, 0];

/// A CRAM writer.
#[derive(Debug)]
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a new CRAM writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let writer = cram::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram as cram;
    /// let writer = cram::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Writes a CRAM file defintion.
    ///
    /// The file ID is set as a blank value (`[0x00; 20]`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram as cram;
    ///
    /// let mut writer = cram::Writer::new(Vec::new());
    /// writer.write_file_definition()?;
    ///
    /// assert_eq!(writer.get_ref(), &[
    ///     // magic number (CRAM)
    ///     0x43, 0x52, 0x41, 0x4d,
    ///     // format (major, minor)
    ///     0x03, 0x00,
    ///     // file ID
    ///     0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ///     0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    /// ]);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_file_definition(&mut self) -> io::Result<()> {
        // magic number
        self.inner.write_all(MAGIC_NUMBER)?;

        // format (major, minor)
        self.inner.write_all(&FILE_DEFINITION_FORMAT)?;

        // File ID is currently blank.
        let file_id = [0; 20];
        self.inner.write_all(&file_id)?;

        Ok(())
    }
}
