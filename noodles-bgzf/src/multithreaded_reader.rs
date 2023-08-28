use std::{
    io::{self, BufRead, Read},
    num::NonZeroUsize,
    thread::{self, JoinHandle},
};

use crossbeam_channel::{Receiver, Sender};

use crate::Block;

type BufferedTx = Sender<io::Result<Block>>;
type BufferedRx = Receiver<io::Result<Block>>;
type InflateTx = Sender<(Vec<u8>, BufferedTx)>;
type InflateRx = Receiver<(Vec<u8>, BufferedTx)>;
type ReadTx = Sender<BufferedRx>;
type ReadRx = Receiver<BufferedRx>;

/// A multithreaded BGZF reader.
///
/// This is a basic multithreaded BGZF reader that uses a thread pool to decompress block data. It
/// differs from a [`super::Reader`] with > 1 worker by placing the inner reader on its own thread
/// to read the raw frames asynchronously.
#[doc(hidden)]
pub struct MultithreadedReader {
    reader_handle: Option<JoinHandle<io::Result<()>>>,
    inflater_handles: Vec<JoinHandle<()>>,
    read_rx: ReadRx,
    position: u64,
    block: Block,
}

impl MultithreadedReader {
    /// Creates a multithreaded BGZF reader.
    pub fn with_worker_count<R>(worker_count: NonZeroUsize, inner: R) -> Self
    where
        R: Read + Send + 'static,
    {
        let (inflate_tx, inflate_rx) = crossbeam_channel::bounded(worker_count.get());
        let (read_tx, read_rx) = crossbeam_channel::bounded(worker_count.get());

        let reader_handle = spawn_reader(inner, inflate_tx, read_tx);
        let inflater_handles = spawn_inflaters(worker_count, inflate_rx);

        Self {
            reader_handle: Some(reader_handle),
            inflater_handles,
            read_rx,
            position: 0,
            block: Block::default(),
        }
    }

    /// Shuts down the reader and inflate workers.
    pub fn finish(&mut self) -> io::Result<()> {
        for handle in self.inflater_handles.drain(..) {
            handle.join().unwrap();
        }

        if let Some(handle) = self.reader_handle.take() {
            handle.join().unwrap()?;
        }

        Ok(())
    }

    fn next_block(&mut self) -> io::Result<Option<Block>> {
        if let Ok(buffered_rx) = self.read_rx.recv() {
            if let Ok(block) = buffered_rx.recv() {
                return block.map(Some);
            }
        }

        Ok(None)
    }

    fn read_block(&mut self) -> io::Result<()> {
        while let Some(mut block) = self.next_block()? {
            block.set_position(self.position);
            self.position += block.size();
            self.block = block;

            if self.block.data().len() > 0 {
                break;
            }
        }

        Ok(())
    }
}

impl Drop for MultithreadedReader {
    fn drop(&mut self) {
        let _ = self.finish();
    }
}

impl Read for MultithreadedReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut src = self.fill_buf()?;
        let amt = src.read(buf)?;
        self.consume(amt);
        Ok(amt)
    }
}

impl BufRead for MultithreadedReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if !self.block.data().has_remaining() {
            self.read_block()?;
        }

        Ok(self.block.data().as_ref())
    }

    fn consume(&mut self, amt: usize) {
        self.block.data_mut().consume(amt);
    }
}

fn spawn_reader<R>(
    mut reader: R,
    inflate_tx: InflateTx,
    read_tx: ReadTx,
) -> JoinHandle<io::Result<()>>
where
    R: Read + Send + 'static,
{
    use super::reader::block::read_frame;

    thread::spawn(move || {
        while let Some(buf) = read_frame(&mut reader)? {
            let (buffered_tx, buffered_rx) = crossbeam_channel::bounded(1);

            if inflate_tx.send((buf, buffered_tx)).is_err() {
                break;
            }

            if read_tx.send(buffered_rx).is_err() {
                break;
            }
        }

        Ok(())
    })
}

fn spawn_inflaters(worker_count: NonZeroUsize, inflate_rx: InflateRx) -> Vec<JoinHandle<()>> {
    use super::reader::block::parse_frame;

    (0..worker_count.get())
        .map(|_| {
            let inflate_rx = inflate_rx.clone();

            thread::spawn(move || {
                while let Ok((src, buffered_tx)) = inflate_rx.recv() {
                    let result = parse_frame(&src);
                    buffered_tx.send(result).unwrap();
                }
            })
        })
        .collect()
}
