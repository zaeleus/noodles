use std::path::Path;

use tokio::{
    fs::File,
    io::{self, AsyncWrite},
};

use super::Writer;
use crate::writer::builder::DEFAULT_LINE_BASE_COUNT;

/// An async FASTA writer builder.
pub struct Builder {
    line_base_count: usize,
}

impl Builder {
    /// Sets the number of bases per line.
    ///
    /// By default, this is set to 80.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::r#async::io::writer::Builder;
    /// let builder = Builder::default().set_line_base_count(100);
    /// ```
    pub fn set_line_base_count(mut self, line_base_count: usize) -> Self {
        self.line_base_count = line_base_count;
        self
    }

    /// Builds an async FASTA writer from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_fasta::r#async::io::writer::Builder;
    /// let writer = Builder::default().build_from_path("in.fa").await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn build_from_path<P>(self, dst: P) -> io::Result<Writer<File>>
    where
        P: AsRef<Path>,
    {
        File::create(dst)
            .await
            .map(|file| self.build_from_writer(file))
    }

    /// Builds an async FASTA writer from an async writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::r#async::io::writer::Builder;
    /// use tokio::io;
    /// let writer = Builder::default().build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> Writer<W>
    where
        W: AsyncWrite + Unpin,
    {
        Writer {
            inner: writer,
            line_base_count: self.line_base_count,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            line_base_count: DEFAULT_LINE_BASE_COUNT,
        }
    }
}
