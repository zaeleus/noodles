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
    Paused(W),
    Running {
        writer_handle: JoinHandle<(W, io::Result<()>)>,
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
    /// let writer = bgzf::io::MultithreadedWriter::new(io::sink());
    /// writer.finish()?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn finish(mut self) -> io::Result<W> {
        self.flush()?;
        self.finish_inner()
    }

    fn finish_inner(&mut self) -> io::Result<W> {
        let state = mem::replace(&mut self.state, State::Done);

        match state {
            State::Paused(inner) => Ok(inner),
            State::Running {
                writer_handle,
                write_tx,
            } => {
                drop(write_tx);
                let (inner, result) = writer_handle.join().unwrap();
                result?;
                Ok(inner)
            }
            State::Done => panic!("invalid state"),
        }
    }

    fn resume(&mut self) {
        if matches!(self.state, State::Running { .. }) {
            return;
        }

        let state = mem::replace(&mut self.state, State::Done);

        let State::Paused(inner) = state else {
            unreachable!();
        };

        let worker_count = self.pool.current_num_threads();
        let (write_tx, write_rx) = crossbeam_channel::bounded(worker_count);
        let writer_handle = spawn_writer(inner, write_rx);

        self.state = State::Running {
            writer_handle,
            write_tx,
        };
    }

    fn remaining(&self) -> usize {
        MAX_BUF_SIZE - self.buf.len()
    }

    fn has_remaining(&self) -> bool {
        self.buf.len() < MAX_BUF_SIZE
    }

    fn send(&mut self) -> io::Result<()> {
        self.resume();

        let State::Running { write_tx, .. } = &self.state else {
            unreachable!();
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
            let _ = self.flush();
            let _ = self.finish_inner();
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

fn spawn_writer<W>(mut writer: W, write_rx: WriteRx) -> JoinHandle<(W, io::Result<()>)>
where
    W: Write + Send + 'static,
{
    thread::spawn(move || {
        let result = write_frames(&mut writer, write_rx);
        (writer, result)
    })
}

fn write_frames<W>(writer: &mut W, write_rx: WriteRx) -> io::Result<()>
where
    W: Write,
{
    use super::writer::{BGZF_EOF, write_frame};

    while let Ok(buffered_rx) = write_rx.recv() {
        if let Ok(result) = buffered_rx.recv() {
            let FrameParts {
                compressed_data,
                crc32,
                uncompressed_size,
            } = result?;

            write_frame(writer, &compressed_data, crc32, uncompressed_size)?;
        }
    }

    writer.write_all(&BGZF_EOF)?;

    Ok(())
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
