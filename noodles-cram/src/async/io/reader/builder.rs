use std::path::Path;

use bytes::BytesMut;
use tokio::{
    fs::File,
    io::{self, AsyncRead},
};

use super::Reader;

/// An async CRAM reader builder.
#[derive(Default)]
pub struct Builder;

impl Builder {
    /// Builds an async CRAM reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_cram::r#async::io::reader::Builder;
    /// use tokio::fs::File;
    /// let _reader = Builder.build_from_path("sample.cram").await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn build_from_path<P>(self, src: P) -> io::Result<Reader<File>>
    where
        P: AsRef<Path>,
    {
        File::open(src)
            .await
            .map(|file| self.build_from_reader(file))
    }

    /// Builds an async CRAM reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::r#async::io::reader::Builder;
    /// use tokio::io;
    /// let _reader = Builder.build_from_reader(io::empty());
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> Reader<R>
    where
        R: AsyncRead + Unpin,
    {
        Reader {
            inner: reader,
            buf: BytesMut::new(),
        }
    }
}
