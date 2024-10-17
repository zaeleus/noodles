use std::path::Path;

use tokio::{
    fs::File,
    io::{self, AsyncWrite},
};

use super::Writer;

/// An async FASTA writer builder.
#[derive(Default)]
pub struct Builder;

impl Builder {
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
        Writer::new(writer)
    }
}
