use std::{
    io::{self, BufRead, Read, Seek, SeekFrom},
    mem,
    num::NonZero,
    thread::{self, JoinHandle},
};

use crossbeam_channel::{Receiver, Sender};

use super::{Block, reader::frame::block_initialize};
use crate::{VirtualPosition, gzi};

type BufferedRx = Receiver<io::Result<Buffer>>;
type ReadTx = Sender<BufferedRx>;
type ReadRx = Receiver<BufferedRx>;
type RecycleTx = Sender<Buffer>;
type RecycleRx = Receiver<Buffer>;

const NON_WORKER_COUNT: usize = 2;

enum State<R> {
    Paused((R, rayon::ThreadPool)),
    Running {
        reader_handle: JoinHandle<Result<(R, rayon::ThreadPool), ReadError<R>>>,
        read_rx: ReadRx,
        recycle_tx: RecycleTx,
    },
    Done,
}

#[derive(Debug, Default)]
struct Buffer {
    buf: Vec<u8>,
    block: Block,
}

/// A multithreaded BGZF reader.
///
/// This is a multithreaded BGZF reader that uses a thread pool to decompress block data. It places
/// the inner reader on its own thread to read raw frames asynchronously.
pub struct MultithreadedReader<R> {
    state: State<R>,
    worker_count: NonZero<usize>,
    position: u64,
    buffer: Buffer,
}

impl<R> MultithreadedReader<R> {
    /// Returns the current position of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// let reader = bgzf::io::MultithreadedReader::new(io::empty());
    /// assert_eq!(reader.position(), 0);
    /// ```
    pub fn position(&self) -> u64 {
        self.position
    }

    /// Returns the current virtual position of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// let reader = bgzf::io::MultithreadedReader::new(io::empty());
    /// assert_eq!(reader.virtual_position(), bgzf::VirtualPosition::MIN);
    /// ```
    pub fn virtual_position(&self) -> VirtualPosition {
        self.buffer.block.virtual_position()
    }

    /// Shuts down the reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// let mut reader = bgzf::io::MultithreadedReader::new(io::empty());
    /// reader.finish()?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn finish(&mut self) -> io::Result<R> {
        let state = mem::replace(&mut self.state, State::Done);

        match state {
            State::Paused((inner, _)) => Ok(inner),
            State::Running {
                reader_handle,
                recycle_tx,
                ..
            } => {
                drop(recycle_tx);

                reader_handle
                    .join()
                    .unwrap()
                    .map(|(inner, _)| inner)
                    .map_err(|e| e.1)
            }
            State::Done => panic!("invalid state"),
        }
    }
}

impl<R> MultithreadedReader<R>
where
    R: Read + Send + 'static,
{
    /// Creates a multithreaded BGZF reader with a worker count of 1.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// let reader = bgzf::io::MultithreadedReader::new(io::empty());
    /// ```
    pub fn new(inner: R) -> Self {
        Self::with_worker_count(NonZero::<usize>::MIN, inner)
    }

    /// Creates a multithreaded BGZF reader with a worker count.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use std::num::NonZero;
    /// use noodles_bgzf as bgzf;
    /// let reader = bgzf::io::MultithreadedReader::with_worker_count(
    ///     NonZero::<usize>::MIN,
    ///     io::empty(),
    /// );
    /// ```
    pub fn with_worker_count(worker_count: NonZero<usize>, inner: R) -> Self {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(worker_count.get())
            .build()
            .unwrap();

        Self {
            state: State::Paused((inner, pool)),
            worker_count,
            position: 0,
            buffer: Buffer::default(),
        }
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// let mut reader = bgzf::io::MultithreadedReader::new(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        self.pause();

        match &mut self.state {
            State::Paused((inner, _)) => inner,
            _ => panic!("invalid state"),
        }
    }

    fn resume(&mut self) {
        if matches!(self.state, State::Running { .. }) {
            return;
        }

        let state = mem::replace(&mut self.state, State::Done);

        let State::Paused((inner, pool)) = state else {
            panic!("invalid state");
        };

        let worker_count = self.worker_count.get();
        let buffer_count = worker_count + NON_WORKER_COUNT;

        let (read_tx, read_rx) = crossbeam_channel::bounded(buffer_count);
        let (recycle_tx, recycle_rx) = crossbeam_channel::bounded(buffer_count);

        for _ in 0..buffer_count {
            recycle_tx.send(Buffer::default()).unwrap();
        }

        let reader_handle = spawn_reader(inner, pool, read_tx, recycle_rx);

        self.state = State::Running {
            reader_handle,
            read_rx,
            recycle_tx,
        };
    }

    fn pause(&mut self) {
        if matches!(self.state, State::Paused(_)) {
            return;
        }

        let state = mem::replace(&mut self.state, State::Done);

        let State::Running {
            reader_handle,
            recycle_tx,
            ..
        } = state
        else {
            panic!("invalid state");
        };

        drop(recycle_tx);

        // Discard read errors.
        let (inner, pool) = match reader_handle.join().unwrap() {
            Ok((inner, pool)) => (inner, pool),
            Err(ReadError((inner, pool), _)) => (inner, pool),
        };

        self.state = State::Paused((inner, pool));
    }

    fn read_block(&mut self) -> io::Result<()> {
        self.resume();

        let State::Running {
            read_rx,
            recycle_tx,
            ..
        } = &self.state
        else {
            panic!("invalid state");
        };

        loop {
            if let Some(mut buffer) = recv_buffer(read_rx)? {
                buffer.block.set_position(self.position);
                self.position += buffer.block.size();

                let prev_buffer = mem::replace(&mut self.buffer, buffer);
                recycle_tx.send(prev_buffer).ok();

                if self.buffer.block.data().len() > 0 {
                    break;
                }
            } else {
                block_initialize(&mut self.buffer.block, 0, 0);
                self.buffer.block.set_position(self.position);
                break;
            }
        }

        Ok(())
    }
}

impl<R> Drop for MultithreadedReader<R> {
    fn drop(&mut self) {
        if !matches!(self.state, State::Done) {
            let _ = self.finish();
        }
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

    fn read_exact(&mut self, buf: &mut [u8]) -> io::Result<()> {
        use super::reader::default_read_exact;

        if let Some(src) = self.buffer.block.data().as_ref().get(..buf.len()) {
            buf.copy_from_slice(src);
            self.consume(src.len());
            Ok(())
        } else {
            default_read_exact(self, buf)
        }
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

impl<R> crate::io::Read for MultithreadedReader<R>
where
    R: Read + Send + 'static,
{
    fn virtual_position(&self) -> VirtualPosition {
        self.buffer.block.virtual_position()
    }
}

impl<R> crate::io::BufRead for MultithreadedReader<R> where R: Read + Send + 'static {}

impl<R> crate::io::Seek for MultithreadedReader<R>
where
    R: Read + Send + Seek + 'static,
{
    fn seek_to_virtual_position(&mut self, pos: VirtualPosition) -> io::Result<VirtualPosition> {
        let (cpos, upos) = pos.into();

        self.get_mut().seek(SeekFrom::Start(cpos))?;
        self.position = cpos;

        self.read_block()?;

        let data = self.buffer.block.data_mut();

        if usize::from(upos) <= data.len() {
            data.set_position(usize::from(upos));
        } else {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }

        Ok(pos)
    }

    fn seek_with_index(&mut self, index: &gzi::Index, pos: SeekFrom) -> io::Result<u64> {
        let SeekFrom::Start(pos) = pos else {
            unimplemented!();
        };

        let virtual_position = index.query(pos)?;
        self.seek_to_virtual_position(virtual_position)?;
        Ok(pos)
    }
}

fn recv_buffer(read_rx: &ReadRx) -> io::Result<Option<Buffer>> {
    if let Ok(buffered_rx) = read_rx.recv()
        && let Ok(buffer) = buffered_rx.recv()
    {
        return buffer.map(Some);
    }

    Ok(None)
}

struct ReadError<R>((R, rayon::ThreadPool), io::Error);

fn spawn_reader<R>(
    mut reader: R,
    pool: rayon::ThreadPool,
    read_tx: ReadTx,
    recycle_rx: RecycleRx,
) -> JoinHandle<Result<(R, rayon::ThreadPool), ReadError<R>>>
where
    R: Read + Send + 'static,
{
    use super::reader::frame::{parse_block, read_frame_into};

    thread::spawn(move || {
        while let Ok(mut buffer) = recycle_rx.recv() {
            match read_frame_into(&mut reader, &mut buffer.buf) {
                Ok(result) if result.is_none() => break,
                Ok(_) => {}
                Err(e) => return Err(ReadError((reader, pool), e)),
            }

            let (buffered_tx, buffered_rx) = crossbeam_channel::bounded(1);

            pool.spawn(move || {
                let result = parse_block(&buffer.buf, &mut buffer.block).map(|_| buffer);
                let _ = buffered_tx.send(result);
            });

            if read_tx.send(buffered_rx).is_err() {
                break;
            }
        }

        Ok((reader, pool))
    })
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;
    use crate::io::Seek;

    #[test]
    fn test_seek_to_virtual_position_with_uncompressed_position_gt_data_len()
    -> Result<(), crate::virtual_position::TryFromU64U16TupleError> {
        #[rustfmt::skip]
        static DATA: &[u8] = &[
            // block 0 (b"noodles")
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x22, 0x00, 0xcb, 0xcb, 0xcf, 0x4f, 0xc9, 0x49, 0x2d, 0x06, 0x00, 0xa1,
            0x58, 0x2a, 0x80, 0x07, 0x00, 0x00, 0x00,
            // EOF block
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        let mut reader = MultithreadedReader::new(Cursor::new(DATA));

        assert!(matches!(
            reader.seek_to_virtual_position(VirtualPosition::try_from((0, 8))?),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_seek_to_virtual_position() -> Result<(), Box<dyn std::error::Error>> {
        use crate::io::Seek;

        #[rustfmt::skip]
        static DATA: &[u8] = &[
            // block 0 (b"noodles")
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x22, 0x00, 0xcb, 0xcb, 0xcf, 0x4f, 0xc9, 0x49, 0x2d, 0x06, 0x00, 0xa1,
            0x58, 0x2a, 0x80, 0x07, 0x00, 0x00, 0x00,
            // EOF block
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        const EOF_VIRTUAL_POSITION: VirtualPosition = VirtualPosition::new(63, 0).unwrap();

        let mut reader = MultithreadedReader::new(Cursor::new(DATA));

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf)?;

        assert_eq!(reader.virtual_position(), EOF_VIRTUAL_POSITION);

        reader.seek_to_virtual_position(const { VirtualPosition::new(0, 3).unwrap() })?;

        buf.clear();
        reader.read_to_end(&mut buf)?;

        assert_eq!(buf, b"dles");
        assert_eq!(reader.virtual_position(), EOF_VIRTUAL_POSITION);

        Ok(())
    }
}
