use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
};

use super::Reader;

const DEFAULT_BUF_CAPACITY: usize = 256 * 1024;

/// A FASTQ reader builder.
///
/// The builder sets a large default buffer capacity (256 KiB) to maximize
/// throughput from the bulk record parsing path.
///
/// # Examples
///
/// ```no_run
/// use noodles_fastq as fastq;
///
/// let mut reader = fastq::io::reader::Builder::default()
///     .build_from_path("reads.fq")?;
///
/// for result in reader.records() {
///     let record = result?;
///     // ...
/// }
/// # Ok::<(), std::io::Error>(())
/// ```
#[derive(Debug)]
pub struct Builder {
    buf_capacity: usize,
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            buf_capacity: DEFAULT_BUF_CAPACITY,
        }
    }
}

impl Builder {
    /// Sets the buffer capacity.
    ///
    /// A larger buffer allows more records to be parsed in bulk from a single
    /// buffer fill, improving throughput. The default is 256 KiB.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_fastq as fastq;
    ///
    /// let mut reader = fastq::io::reader::Builder::default()
    ///     .buf_capacity(512 * 1024)
    ///     .build_from_path("reads.fq")?;
    /// # Ok::<(), std::io::Error>(())
    /// ```
    pub fn buf_capacity(mut self, buf_capacity: usize) -> Self {
        self.buf_capacity = buf_capacity;
        self
    }

    /// Builds a FASTQ reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_fastq as fastq;
    ///
    /// let mut reader = fastq::io::reader::Builder::default()
    ///     .build_from_path("reads.fq")?;
    /// # Ok::<(), std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, src: P) -> io::Result<Reader<BufReader<File>>>
    where
        P: AsRef<Path>,
    {
        let file = File::open(src)?;
        let buf_reader = BufReader::with_capacity(self.buf_capacity, file);
        Ok(Reader::new(buf_reader))
    }

    /// Builds a FASTQ reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    ///
    /// let data = b"@r0\nACGT\n+\nNDLS\n";
    /// let mut reader = fastq::io::reader::Builder::default()
    ///     .build_from_reader(&data[..]);
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> Reader<R>
    where
        R: BufRead,
    {
        Reader::new(reader)
    }
}
