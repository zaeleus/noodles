//! Multithreaded BGZF writer.

mod builder;

use std::{
    io::{self, Write},
    mem,
    num::NonZero,
    sync::{
        Arc,
        atomic::{AtomicU64, Ordering},
    },
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
type WriteTx = Sender<(BufferedRx, u64)>; // Now includes block_number
type WriteRx = Receiver<(BufferedRx, u64)>;
type BlockInfoTx = Sender<BlockInfo>;

/// Block completion info sent from writer thread.
///
/// This enables building indexes while using multi-threaded compression:
/// cache index entries with their block number, then resolve positions
/// when the block completion notification is received.
#[derive(Debug, Clone, Copy)]
pub struct BlockInfo {
    /// The block number (0-indexed, assigned when block is sent for compression).
    pub block_number: u64,
    /// The compressed file position where this block starts.
    pub compressed_start: u64,
    /// The size of the compressed block (including header and trailer).
    pub compressed_size: usize,
    /// The size of the uncompressed data in this block.
    pub uncompressed_size: usize,
}

/// Receiver for block completion notifications.
pub type BlockInfoRx = Receiver<BlockInfo>;

enum State<W> {
    Running {
        writer_handle: JoinHandle<io::Result<W>>,
        deflater_handles: Vec<JoinHandle<()>>,
        write_tx: WriteTx,
        deflate_tx: DeflateTx,
        block_info_rx: BlockInfoRx,
    },
    Done,
}

/// A multithreaded BGZF writer.
///
/// This is much more basic than [`super::Writer`] but uses a thread pool to compress block data.
///
/// # Position Tracking
///
/// This writer supports position tracking for building BAM indexes during multi-threaded
/// compression. Use [`block_info_receiver`](Self::block_info_receiver) to receive
/// notifications when blocks are written to the output, then resolve cached index entries
/// with their final compressed positions.
pub struct MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    state: State<W>,
    buf: BytesMut,
    // Block tracking for indexing
    current_block_number: u64,
    position: Arc<AtomicU64>,
    blocks_written: Arc<AtomicU64>,
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

        let state = mem::replace(&mut self.state, State::Done);

        match state {
            State::Running {
                writer_handle,
                mut deflater_handles,
                write_tx,
                deflate_tx,
                block_info_rx: _,
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

        // Include block number with the send
        write_tx
            .send((buffered_rx, self.current_block_number))
            .unwrap();
        self.current_block_number += 1;

        let src = self.buf.split().freeze();
        let message = (src, buffered_tx);
        deflate_tx.send(message).unwrap();

        Ok(())
    }

    /// Returns a receiver for block completion notifications.
    ///
    /// Use this for building indexes: cache index entries with their block number,
    /// then resolve positions when the block completion notification is received.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io::{self, Write};
    /// use noodles_bgzf as bgzf;
    ///
    /// let mut writer = bgzf::io::MultithreadedWriter::new(io::sink());
    ///
    /// // Get the receiver before writing
    /// let rx = writer.block_info_receiver().unwrap().clone();
    ///
    /// // Write data
    /// writer.write_all(b"hello")?;
    /// writer.flush()?;
    ///
    /// // Check for block completions (non-blocking)
    /// while let Ok(info) = rx.try_recv() {
    ///     println!("Block {} written at position {}", info.block_number, info.compressed_start);
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn block_info_receiver(&self) -> Option<&BlockInfoRx> {
        match &self.state {
            State::Running { block_info_rx, .. } => Some(block_info_rx),
            State::Done => None,
        }
    }

    /// Returns the current block number (the next block to be written).
    ///
    /// This is incremented each time data is flushed (a new block is sent for compression).
    pub fn current_block_number(&self) -> u64 {
        self.current_block_number
    }

    /// Returns the number of blocks fully written to the output.
    ///
    /// This lags behind `current_block_number` while blocks are being compressed
    /// and written.
    pub fn blocks_written(&self) -> u64 {
        self.blocks_written.load(Ordering::Acquire)
    }

    /// Returns the current compressed file position.
    ///
    /// This is the total number of bytes written to the underlying writer.
    pub fn position(&self) -> u64 {
        self.position.load(Ordering::Acquire)
    }

    /// Returns the current offset in the staging buffer (uncompressed bytes).
    ///
    /// This is the number of bytes written since the last flush.
    pub fn buffer_offset(&self) -> usize {
        self.buf.len()
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

fn spawn_writer<W>(
    mut writer: W,
    write_rx: WriteRx,
    position: Arc<AtomicU64>,
    blocks_written: Arc<AtomicU64>,
    block_info_tx: BlockInfoTx,
) -> JoinHandle<io::Result<W>>
where
    W: Write + Send + 'static,
{
    use super::writer::{BGZF_EOF, write_frame};

    thread::spawn(move || {
        while let Ok((buffered_rx, block_number)) = write_rx.recv() {
            if let Ok(result) = buffered_rx.recv() {
                let (compressed_data, crc32, uncompressed_size) = result?;

                // Capture position before writing
                let compressed_start = position.load(Ordering::Acquire);

                // Write the block
                let block_size =
                    write_frame(&mut writer, &compressed_data, crc32, uncompressed_size)?;

                // Update position after writing
                position.fetch_add(block_size as u64, Ordering::Release);
                blocks_written.fetch_add(1, Ordering::Release);

                // Send block completion info
                let _ = block_info_tx.send(BlockInfo {
                    block_number,
                    compressed_start,
                    compressed_size: block_size,
                    uncompressed_size,
                });
            }
        }

        writer.write_all(&BGZF_EOF)?;
        position.fetch_add(BGZF_EOF.len() as u64, Ordering::Release);

        Ok(writer)
    })
}

fn spawn_deflaters<L>(
    compression_level: L,
    worker_count: NonZero<usize>,
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
