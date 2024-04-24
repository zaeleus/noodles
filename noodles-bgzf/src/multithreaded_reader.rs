use std::{
    io::{self, BufRead, Read},
    mem,
    num::NonZeroUsize,
    thread::{self, JoinHandle},
};

use crossbeam_channel::{Receiver, Sender};

use crate::{Block, VirtualPosition};

type BufferedTx = Sender<io::Result<Buffer>>;
type BufferedRx = Receiver<io::Result<Buffer>>;
type InflateTx = Sender<(Buffer, BufferedTx)>;
type InflateRx = Receiver<(Buffer, BufferedTx)>;
type ReadTx = Sender<BufferedRx>;
type ReadRx = Receiver<BufferedRx>;
type RecycleTx = Sender<Buffer>;
type RecycleRx = Receiver<Buffer>;

enum State<R> {
    Paused(R),
    Running {
        reader_handle: JoinHandle<Result<R, ReadError<R>>>,
        inflater_handles: Vec<JoinHandle<()>>,
        read_rx: ReadRx,
        recycle_tx: RecycleTx,
    },
}

#[derive(Debug, Default)]
struct Buffer {
    buf: Vec<u8>,
    block: Block,
}

/// A multithreaded BGZF reader.
///
/// This is a basic multithreaded BGZF reader that uses a thread pool to decompress block data. It
/// differs from a [`super::Reader`] with > 1 worker by placing the inner reader on its own thread
/// to read the raw frames asynchronously.
pub struct MultithreadedReader<R> {
    state: Option<State<R>>,
    worker_count: NonZeroUsize,
    position: u64,
    buffer: Buffer,
}

impl<R> MultithreadedReader<R> {
    /// Returns the current position of the stream.
    pub fn position(&self) -> u64 {
        self.position
    }

    /// Returns the current virtual position of the stream.
    pub fn virtual_position(&self) -> VirtualPosition {
        self.buffer.block.virtual_position()
    }

    /// Shuts down the reader and inflate workers.
    pub fn finish(&mut self) -> io::Result<R> {
        let State::Running {
            reader_handle,
            mut inflater_handles,
            recycle_tx,
            ..
        } = self.state.take().unwrap()
        else {
            panic!("invalid state");
        };

        drop(recycle_tx);

        for handle in inflater_handles.drain(..) {
            handle.join().unwrap();
        }

        reader_handle.join().unwrap().map_err(|e| e.1)
    }

    fn recv_buffer(&self) -> io::Result<Option<Buffer>> {
        let State::Running { read_rx, .. } = self.state.as_ref().unwrap() else {
            panic!("invalid state");
        };

        if let Ok(buffered_rx) = read_rx.recv() {
            if let Ok(buffer) = buffered_rx.recv() {
                return buffer.map(Some);
            }
        }

        Ok(None)
    }
}

impl<R> MultithreadedReader<R>
where
    R: Read + Send + 'static,
{
    /// Creates a multithreaded BGZF reader.
    pub fn with_worker_count(worker_count: NonZeroUsize, inner: R) -> Self {
        Self {
            state: Some(State::Paused(inner)),
            worker_count,
            position: 0,
            buffer: Buffer::default(),
        }
    }

    fn resume(&mut self) {
        if matches!(self.state, Some(State::Running { .. })) {
            return;
        }

        let Some(State::Paused(inner)) = self.state.take() else {
            panic!("invalid state");
        };

        let worker_count = self.worker_count.get();

        let (inflate_tx, inflate_rx) = crossbeam_channel::bounded(worker_count);
        let (read_tx, read_rx) = crossbeam_channel::bounded(worker_count);
        let (recycle_tx, recycle_rx) = crossbeam_channel::bounded(worker_count);

        for _ in 0..worker_count {
            recycle_tx.send(Buffer::default()).unwrap();
        }

        let reader_handle = spawn_reader(inner, inflate_tx, read_tx, recycle_rx);
        let inflater_handles = spawn_inflaters(self.worker_count, inflate_rx);

        self.state = Some(State::Running {
            reader_handle,
            inflater_handles,
            read_rx,
            recycle_tx,
        });
    }

    fn read_block(&mut self) -> io::Result<()> {
        self.resume();

        let State::Running { recycle_tx, .. } = self.state.as_ref().unwrap() else {
            panic!("invalid state");
        };

        while let Some(mut buffer) = self.recv_buffer()? {
            buffer.block.set_position(self.position);
            self.position += buffer.block.size();

            let prev_buffer = mem::replace(&mut self.buffer, buffer);
            recycle_tx.send(prev_buffer).ok();

            if self.buffer.block.data().len() > 0 {
                break;
            }
        }

        Ok(())
    }
}

impl<R> Drop for MultithreadedReader<R> {
    fn drop(&mut self) {
        let _ = self.finish();
    }
}

impl<R> Read for MultithreadedReader<R>
where
    R: Read + Send + 'static,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut src = self.fill_buf()?;
        let amt = src.read(buf)?;
        self.consume(amt);
        Ok(amt)
    }
}

impl<R> BufRead for MultithreadedReader<R>
where
    R: Read + Send + 'static,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if !self.buffer.block.data().has_remaining() {
            self.read_block()?;
        }

        Ok(self.buffer.block.data().as_ref())
    }

    fn consume(&mut self, amt: usize) {
        self.buffer.block.data_mut().consume(amt);
    }
}

struct ReadError<R>(R, io::Error);

fn spawn_reader<R>(
    mut reader: R,
    inflate_tx: InflateTx,
    read_tx: ReadTx,
    recycle_rx: RecycleRx,
) -> JoinHandle<Result<R, ReadError<R>>>
where
    R: Read + Send + 'static,
{
    use super::reader::block::read_frame_into;

    thread::spawn(move || {
        while let Ok(mut buffer) = recycle_rx.recv() {
            match read_frame_into(&mut reader, &mut buffer.buf) {
                Ok(result) if result.is_none() => break,
                Ok(_) => {}
                Err(e) => return Err(ReadError(reader, e)),
            }

            let (buffered_tx, buffered_rx) = crossbeam_channel::bounded(1);

            inflate_tx.send((buffer, buffered_tx)).unwrap();
            read_tx.send(buffered_rx).unwrap();
        }

        Ok(reader)
    })
}

fn spawn_inflaters(worker_count: NonZeroUsize, inflate_rx: InflateRx) -> Vec<JoinHandle<()>> {
    use super::reader::block::parse_block;

    (0..worker_count.get())
        .map(|_| {
            let inflate_rx = inflate_rx.clone();

            thread::spawn(move || {
                while let Ok((mut buffer, buffered_tx)) = inflate_rx.recv() {
                    let result = parse_block(&buffer.buf, &mut buffer.block).map(|_| buffer);
                    buffered_tx.send(result).unwrap();
                }
            })
        })
        .collect()
}
