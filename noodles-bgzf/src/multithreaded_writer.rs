//! Multithreaded BGZF writer.

mod builder;

use std::{
    io::{self, Write},
    mem,
    num::NonZeroUsize,
    thread::{self, JoinHandle},
};

use byteorder::{LittleEndian, WriteBytesExt};
use bytes::{Bytes, BytesMut};
use crossbeam_channel::{Receiver, Sender};

pub use self::builder::Builder;
use super::{gz, writer::CompressionLevelImpl};

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
    /// let writer = bgzf::MultithreadedWriter::new(io::sink());
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
    /// let writer = bgzf::MultithreadedWriter::with_worker_count(NonZeroUsize::MIN, io::sink());
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
    /// let mut writer = bgzf::MultithreadedWriter::new(io::sink());
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
                let (compressed_data, crc32, uncompressed_len) = result?;
                write_frame(&mut writer, &compressed_data, crc32, uncompressed_len)?;
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

fn compress(src: &[u8], compression_level: CompressionLevelImpl) -> io::Result<FrameParts> {
    use super::deflate;
    let (cdata, crc32, _) = deflate::encode(src, compression_level)?;
    Ok((cdata, crc32, src.len()))
}

fn write_frame<W>(
    writer: &mut W,
    compressed_data: &[u8],
    crc32: u32,
    uncompressed_len: usize,
) -> io::Result<()>
where
    W: Write,
{
    use super::BGZF_HEADER_SIZE;

    let block_size = BGZF_HEADER_SIZE + compressed_data.len() + gz::TRAILER_SIZE;
    write_header(writer, block_size)?;

    writer.write_all(compressed_data)?;

    write_trailer(writer, crc32, uncompressed_len)?;

    Ok(())
}

fn write_header<W>(writer: &mut W, block_size: usize) -> io::Result<()>
where
    W: Write,
{
    const BGZF_FLG: u8 = 0x04; // FEXTRA
    const BGZF_XFL: u8 = 0x00; // none
    const BGZF_XLEN: u16 = 6;

    const BGZF_SI1: u8 = b'B';
    const BGZF_SI2: u8 = b'C';
    const BGZF_SLEN: u16 = 2;

    writer.write_all(&gz::MAGIC_NUMBER)?;
    writer.write_u8(gz::CompressionMethod::Deflate as u8)?;
    writer.write_u8(BGZF_FLG)?;
    writer.write_u32::<LittleEndian>(gz::MTIME_NONE)?;
    writer.write_u8(BGZF_XFL)?;
    writer.write_u8(gz::OperatingSystem::Unknown as u8)?;
    writer.write_u16::<LittleEndian>(BGZF_XLEN)?;

    writer.write_u8(BGZF_SI1)?;
    writer.write_u8(BGZF_SI2)?;
    writer.write_u16::<LittleEndian>(BGZF_SLEN)?;

    let bsize = u16::try_from(block_size - 1)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u16::<LittleEndian>(bsize)?;

    Ok(())
}

fn write_trailer<W>(writer: &mut W, crc32: u32, uncompressed_len: usize) -> io::Result<()>
where
    W: Write,
{
    writer.write_u32::<LittleEndian>(crc32)?;

    let r#isize = u32::try_from(uncompressed_len)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32::<LittleEndian>(r#isize)?;

    Ok(())
}
