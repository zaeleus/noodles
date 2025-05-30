//! Multithreaded BGZF writer.

mod builder;

use std::{
    io::{self, Write},
    mem,
    num::NonZeroUsize,
    thread::{self, JoinHandle},
};

use bytes::{Bytes, BytesMut};
use crossbeam_channel::{Receiver, Sender};

pub use self::builder::Builder;
use super::writer::{CompressionLevelImpl, MAX_BUF_SIZE};

type FrameParts = (Vec<u8>, u32, usize);
type BufferedTx = Sender<io::Result<FrameParts>>;
type BufferedRx = Receiver<io::Result<FrameParts>>;
type DeflateTx = Sender<(Bytes, BufferedTx)>;
type DeflateRx = Receiver<(Bytes, BufferedTx)>;
type WriteTx = Sender<BufferedRx>;
type WriteRx = Receiver<BufferedRx>;

enum State<W> {
    Running {
        writer_handle: JoinHandle<io::Result<W>>,
        deflater_handles: Vec<JoinHandle<()>>,
        write_tx: WriteTx,
        deflate_tx: DeflateTx,
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
    state: State<W>,
    buf: BytesMut,
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
        Builder::default().build_from_writer(inner)
    }

    /// Creates a multithreaded BGZF writer with a worker count.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use std::num::NonZeroUsize;
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::io::MultithreadedWriter::with_worker_count(NonZeroUsize::MIN, io::sink());
    /// ```
    pub fn with_worker_count(worker_count: NonZeroUsize, inner: W) -> Self {
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

        let state = mem::replace(&mut self.state, State::Done);

        match state {
            State::Running {
                writer_handle,
                mut deflater_handles,
                write_tx,
                deflate_tx,
            } => {
                drop(deflate_tx);

                for handle in deflater_handles.drain(..) {
                    handle.join().unwrap();
                }

                drop(write_tx);

                writer_handle.join().unwrap()
            }
            State::Done => panic!("invalid state"),
        }
    }

    fn remaining(&self) -> usize {
        MAX_BUF_SIZE - self.buf.len()
    }

    fn has_remaining(&self) -> bool {
        self.buf.len() < MAX_BUF_SIZE
    }

    fn send(&mut self) -> io::Result<()> {
        let State::Running {
            write_tx,
            deflate_tx,
            ..
        } = &self.state
        else {
            panic!("invalid state");
        };

        let (buffered_tx, buffered_rx) = crossbeam_channel::bounded(1);

        write_tx.send(buffered_rx).unwrap();

        let src = self.buf.split().freeze();
        let message = (src, buffered_tx);
        deflate_tx.send(message).unwrap();

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

fn spawn_writer<W>(mut writer: W, write_rx: WriteRx) -> JoinHandle<io::Result<W>>
where
    W: Write + Send + 'static,
{
    use super::writer::{BGZF_EOF, write_frame};

    thread::spawn(move || {
        while let Ok(buffered_rx) = write_rx.recv() {
            if let Ok(result) = buffered_rx.recv() {
                let (compressed_data, crc32, uncompressed_size) = result?;
                write_frame(&mut writer, &compressed_data, crc32, uncompressed_size)?;
            }
        }

        writer.write_all(&BGZF_EOF)?;

        Ok(writer)
    })
}

fn spawn_deflaters<L>(
    compression_level: L,
    worker_count: NonZeroUsize,
    deflate_rx: DeflateRx,
) -> Vec<JoinHandle<()>>
where
    L: Into<CompressionLevelImpl>,
{
    let compression_level = compression_level.into();

    (0..worker_count.get())
        .map(|_| {
            let deflate_rx = deflate_rx.clone();

            thread::spawn(move || {
                while let Ok((src, buffered_tx)) = deflate_rx.recv() {
                    let result = compress(&src, compression_level);
                    buffered_tx.send(result).ok();
                }
            })
        })
        .collect()
}

fn compress(src: &[u8], compression_level: CompressionLevelImpl) -> io::Result<FrameParts> {
    use crate::deflate;

    let mut dst = Vec::new();
    let crc32 = deflate::encode(src, compression_level, &mut dst)?;
    Ok((dst, crc32, src.len()))
}
