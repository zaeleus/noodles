//! Multithreaded BGZF writer.

mod builder;

use std::{
    io::{self, Write},
    mem,
    num::NonZero,
    thread::{self, JoinHandle},
};

use bytes::BytesMut;
use crossbeam_channel::{Receiver, Sender};

pub use self::builder::Builder;
use super::writer::{CompressionLevelImpl, MAX_BUF_SIZE};

type BufferedRx = Receiver<io::Result<FrameParts>>;
type WriteTx = Sender<BufferedRx>;
type WriteRx = Receiver<BufferedRx>;

struct FrameParts {
    compressed_data: Vec<u8>,
    crc32: u32,
    uncompressed_size: usize,
}

enum State<W> {
    Running {
        writer_handle: JoinHandle<io::Result<W>>,
        write_tx: WriteTx,
    },
    Done,
}

/// A multithreaded BGZF writer.
///
/// This is much more basic than [`super::Writer`] but uses a thread pool to compress block data.
pub struct MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    pool: rayon::ThreadPool,
    state: State<W>,
    buf: BytesMut,
    compression_level: CompressionLevelImpl,
}

impl<W> MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    /// Creates a multithreaded BGZF writer with a default worker count.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::io::MultithreadedWriter::new(io::sink());
    /// ```
    pub fn new(inner: W) -> Self {
        Self::with_worker_count(NonZero::<usize>::MIN, inner)
    }

    /// Creates a multithreaded BGZF writer with a worker count.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use std::num::NonZero;
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::io::MultithreadedWriter::with_worker_count(
    ///     NonZero::<usize>::MIN,
    ///     io::sink(),
    /// );
    /// ```
    pub fn with_worker_count(worker_count: NonZero<usize>, inner: W) -> Self {
        Builder::default()
            .set_worker_count(worker_count)
            .build_from_writer(inner)
    }

    /// Finishes the output stream by flushing any remaining buffers.
    ///
    /// This shuts down the writer and deflater workers and appends the final BGZF EOF block.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// let mut writer = bgzf::io::MultithreadedWriter::new(io::sink());
    /// writer.finish()?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn finish(&mut self) -> io::Result<W> {
        self.flush()?;
        self.finish_inner()
    }

    fn finish_inner(&mut self) -> io::Result<W> {
        let state = mem::replace(&mut self.state, State::Done);

        match state {
            State::Running {
                writer_handle,
                write_tx,
            } => {
                drop(write_tx);
                writer_handle.join().unwrap()
            }
            State::Done => Err(finished_error()),
        }
    }

    fn remaining(&self) -> usize {
        MAX_BUF_SIZE - self.buf.len()
    }

    fn has_remaining(&self) -> bool {
        self.buf.len() < MAX_BUF_SIZE
    }

    fn send(&mut self) -> io::Result<()> {
        let State::Running { write_tx, .. } = &self.state else {
            return Err(finished_error());
        };

        let (buffered_tx, buffered_rx) = crossbeam_channel::bounded(1);

        if write_tx.send(buffered_rx).is_err() {
            return self.finish_inner().map(|_| ());
        }

        let src = self.buf.split().freeze();
        let compression_level = self.compression_level;

        self.pool.spawn(move || {
            let result = compress(&src, compression_level);
            buffered_tx.send(result).ok();
        });

        Ok(())
    }
}

impl<W> Drop for MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    fn drop(&mut self) {
        if !matches!(self.state, State::Done) {
            let _ = self.finish();
        }
    }
}

impl<W> Write for MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let amt = self.remaining().min(buf.len());
        self.buf.extend_from_slice(&buf[..amt]);

        if !self.has_remaining() {
            self.flush()?;
        }

        Ok(amt)
    }

    fn flush(&mut self) -> io::Result<()> {
        if self.buf.is_empty() {
            Ok(())
        } else {
            self.send()
        }
    }
}

// The writer is finished after `finish` or after a write error, which is returned once.
fn finished_error() -> io::Error {
    io::Error::other("multithreaded writer is finished")
}

fn spawn_writer<W>(mut writer: W, write_rx: WriteRx) -> JoinHandle<io::Result<W>>
where
    W: Write + Send + 'static,
{
    use super::writer::{BGZF_EOF, write_frame};

    thread::spawn(move || {
        while let Ok(buffered_rx) = write_rx.recv() {
            if let Ok(result) = buffered_rx.recv() {
                let FrameParts {
                    compressed_data,
                    crc32,
                    uncompressed_size,
                } = result?;

                write_frame(&mut writer, &compressed_data, crc32, uncompressed_size)?;
            }
        }

        writer.write_all(&BGZF_EOF)?;

        Ok(writer)
    })
}

fn compress(src: &[u8], compression_level: CompressionLevelImpl) -> io::Result<FrameParts> {
    use crate::deflate;

    let mut dst = Vec::new();
    let crc32 = deflate::encode(src, compression_level, &mut dst)?;

    Ok(FrameParts {
        compressed_data: dst,
        crc32,
        uncompressed_size: src.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    struct FailingWriter;

    impl Write for FailingWriter {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::from(io::ErrorKind::StorageFull))
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn test_write_and_finish_after_write_error() {
        let mut writer = MultithreadedWriter::new(FailingWriter);

        // The write error is returned asynchronously, by a later write or by finish.
        let mut result = Ok(());

        for _ in 0..64 {
            result = writer.write_all(&[0; MAX_BUF_SIZE]);

            if result.is_err() {
                break;
            }
        }

        if result.is_ok() {
            result = writer.finish().map(|_| ());
        }

        assert!(matches!(result, Err(e) if e.kind() == io::ErrorKind::StorageFull));

        assert!(writer.write_all(&[0; MAX_BUF_SIZE]).is_err());
        assert!(writer.finish().is_err());
    }

    #[test]
    fn test_finish_twice() -> io::Result<()> {
        let mut writer = MultithreadedWriter::new(io::sink());
        writer.finish()?;
        assert!(writer.finish().is_err());
        Ok(())
    }
}
