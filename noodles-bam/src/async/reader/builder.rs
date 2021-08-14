use noodles_bgzf as bgzf;
use tokio::io::AsyncRead;

use super::Reader;

/// An async BAM reader builder.
pub struct Builder<R> {
    inner: R,
    worker_count: Option<usize>,
}

impl<R> Builder<R>
where
    R: AsyncRead + Unpin,
{
    pub(crate) fn new(inner: R) -> Self {
        Self {
            inner,
            worker_count: None,
        }
    }

    /// Sets a worker count.
    ///
    /// By default, the worker count is set to the number of available logical CPUs.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let data = [];
    /// let builder = bam::AsyncReader::builder(&data[..]).set_worker_count(8);
    /// ```
    pub fn set_worker_count(mut self, worker_count: usize) -> Self {
        self.worker_count = Some(worker_count);
        self
    }

    /// Builds an async BAM reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let data = [];
    /// let reader = bam::AsyncReader::builder(&data[..]).build();
    /// ```
    pub fn build(self) -> Reader<R> {
        let mut builder = bgzf::AsyncReader::builder(self.inner);

        if let Some(worker_count) = self.worker_count {
            builder = builder.set_worker_count(worker_count);
        }

        Reader {
            inner: builder.build(),
        }
    }
}
