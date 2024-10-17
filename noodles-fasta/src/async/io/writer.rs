//! Async FASTA writer.

mod builder;
mod record;

use tokio::io::{self, AsyncWrite};

pub use self::builder::Builder;
use self::record::write_record;
use crate::Record;

/// An async FASTA writer.
pub struct Writer<W> {
    inner: W,
    line_base_count: usize,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// use tokio::io;
    /// let writer = fasta::r#async::io::Writer::new(io::sink());
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
    /// use noodles_fasta as fasta;
    /// use tokio::io;
    /// let mut writer = fasta::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Unwraps and returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// use tokio::io;
    /// let writer = fasta::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }
}

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Creates a FASTA writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// use tokio::io;
    /// let writer = fasta::r#async::io::Writer::new(io::sink());
    /// ```
    pub fn new(inner: W) -> Self {
        Builder::default().build_from_writer(inner)
    }

    /// Writes a FASTA record.
    ///
    /// By default, sequence lines are hard wrapped at 80 bases. This can be changed by using
    /// [`Builder::set_line_base_count`] when creating the writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_fasta::{self as fasta, record::{Definition, Sequence}};
    /// use tokio::io;
    ///
    /// let mut writer = fasta::r#async::io::Writer::new(io::sink());
    ///
    /// let definition = Definition::new("sq0", None);
    /// let sequence = Sequence::from(b"ACGT".to_vec());
    /// let record = fasta::Record::new(definition, sequence);
    ///
    /// writer.write_record(&record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_record(&mut self, record: &Record) -> io::Result<()> {
        write_record(&mut self.inner, record, self.line_base_count).await
    }
}
