//! Multithreaded BGZF writer.

mod builder;

use std::{
    io::{self, Write},
    num::NonZeroUsize,
    thread::{self, JoinHandle},
};

use bytes::{BufMut, Bytes, BytesMut};
use crossbeam_channel::{Receiver, Sender};

pub use self::builder::Builder;
use super::{
    gz,
    writer::{CompressionLevel, CompressionLevelImpl},
};

type BufferedTx = Sender<io::Result<Vec<u8>>>;
type BufferedRx = Receiver<io::Result<Vec<u8>>>;
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
}

/// A multithreaded BGZF writer.
///
/// This is much more basic than [`super::Writer`] but uses a thread pool to compress block data.
pub struct MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    state: Option<State<W>>,
    buf: BytesMut,
}

impl<W> MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    fn with_compression_level_and_worker_count(
        compression_level: CompressionLevel,
        worker_count: NonZeroUsize,
        inner: W,
    ) -> Self {
        let (write_tx, write_rx) = crossbeam_channel::bounded(worker_count.get());
        let (deflate_tx, deflate_rx) = crossbeam_channel::bounded(worker_count.get());

        let writer_handle = spawn_writer(inner, write_rx);
        let deflater_handles = spawn_deflaters(compression_level, worker_count, deflate_rx);

        Self {
            state: Some(State::Running {
                writer_handle,
                deflater_handles,
                write_tx,
                deflate_tx,
            }),
            buf: BytesMut::new(),
        }
    }

    /// Creates a multithreaded BGZF writer with a default worker count.
    pub fn new(inner: W) -> Self {
        Self::with_worker_count(NonZeroUsize::MIN, inner)
    }

    /// Creates a multithreaded BGZF writer with a worker count.
    pub fn with_worker_count(worker_count: NonZeroUsize, inner: W) -> Self {
        Self::with_compression_level_and_worker_count(
            CompressionLevel::default(),
            worker_count,
            inner,
        )
    }

    /// Finishes the output stream by flushing any remaining buffers.
    ///
    /// This shuts down the writer and deflater workers and appends the final BGZF EOF block.
    pub fn finish(&mut self) -> io::Result<W> {
        self.flush()?;

        match self.state.take().unwrap() {
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
        }
    }

    fn send(&mut self) -> io::Result<()> {
        let State::Running {
            write_tx,
            deflate_tx,
            ..
        } = self.state.as_ref().unwrap();

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
        if self.state.is_some() {
            let _ = self.finish();
        }
    }
}

impl<W> Write for MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        use std::cmp;

        use super::writer::MAX_BUF_SIZE;

        let amt = cmp::min(MAX_BUF_SIZE - self.buf.len(), buf.len());
        self.buf.extend_from_slice(&buf[..amt]);

        if self.buf.len() >= MAX_BUF_SIZE {
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
    use super::writer::BGZF_EOF;

    thread::spawn(move || {
        while let Ok(buffered_rx) = write_rx.recv() {
            if let Ok(result) = buffered_rx.recv() {
                let buf = result?;
                writer.write_all(&buf[..])?;
            }
        }

        writer.write_all(BGZF_EOF)?;

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

fn compress(src: &[u8], compression_level: CompressionLevelImpl) -> io::Result<Vec<u8>> {
    use super::{writer::deflate_data, BGZF_HEADER_SIZE};

    let mut dst = Vec::new();

    let (cdata, crc32, _) = deflate_data(src, compression_level)?;

    let block_size = BGZF_HEADER_SIZE + cdata.len() + gz::TRAILER_SIZE;
    put_header(&mut dst, block_size)?;

    dst.extend(cdata);

    put_trailer(&mut dst, crc32, src.len())?;

    Ok(dst)
}

fn put_header<B>(dst: &mut B, block_size: usize) -> io::Result<()>
where
    B: BufMut,
{
    const BGZF_FLG: u8 = 0x04; // FEXTRA
    const BGZF_XFL: u8 = 0x00; // none
    const BGZF_XLEN: u16 = 6;

    const BGZF_SI1: u8 = b'B';
    const BGZF_SI2: u8 = b'C';
    const BGZF_SLEN: u16 = 2;

    dst.put_slice(&gz::MAGIC_NUMBER);
    dst.put_u8(gz::CompressionMethod::Deflate as u8);
    dst.put_u8(BGZF_FLG);
    dst.put_u32_le(gz::MTIME_NONE);
    dst.put_u8(BGZF_XFL);
    dst.put_u8(gz::OperatingSystem::Unknown as u8);
    dst.put_u16_le(BGZF_XLEN);

    dst.put_u8(BGZF_SI1);
    dst.put_u8(BGZF_SI2);
    dst.put_u16_le(BGZF_SLEN);

    let bsize = u16::try_from(block_size - 1)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.put_u16_le(bsize);

    Ok(())
}

fn put_trailer<B>(dst: &mut B, crc32: u32, uncompressed_len: usize) -> io::Result<()>
where
    B: BufMut,
{
    dst.put_u32_le(crc32);

    let r#isize = u32::try_from(uncompressed_len)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.put_u32_le(r#isize);

    Ok(())
}
