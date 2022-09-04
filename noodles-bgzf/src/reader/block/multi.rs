use std::{
    collections::VecDeque,
    io::{self, Read},
    num::NonZeroUsize,
    thread::{self, JoinHandle},
};

use crossbeam_channel::{Receiver, Sender};

use crate::Block;

type BufferedTx = Sender<io::Result<Block>>;
type BufferedRx = Receiver<io::Result<Block>>;
type InflaterTx = Sender<(Vec<u8>, BufferedTx)>;
type InflaterRx = Receiver<(Vec<u8>, BufferedTx)>;

pub struct Reader<R> {
    inner: Option<R>,
    inflater_tx: Option<InflaterTx>,
    inflater_handles: Vec<JoinHandle<()>>,
    queue: VecDeque<BufferedRx>,
    is_eof: bool,
}

impl<R> Reader<R> {
    fn shutdown(&mut self) -> io::Result<()> {
        self.inflater_tx.take();

        for handle in self.inflater_handles.drain(..) {
            handle.join().unwrap();
        }

        Ok(())
    }
}

impl<R> Drop for Reader<R> {
    fn drop(&mut self) {
        let _ = self.shutdown();
    }
}

impl<R> Reader<R>
where
    R: Read,
{
    pub(crate) fn with_worker_count(worker_count: NonZeroUsize, inner: R) -> Self {
        let worker_count = worker_count.get();

        let (inflater_tx, inflater_rx) = crossbeam_channel::bounded(worker_count);
        let inflater_handles = spawn_inflaters(worker_count, inflater_rx);

        Self {
            inner: Some(inner),
            inflater_tx: Some(inflater_tx),
            inflater_handles,
            queue: VecDeque::with_capacity(worker_count),
            is_eof: false,
        }
    }

    pub fn get_ref(&self) -> &R {
        self.inner.as_ref().unwrap()
    }

    pub fn get_mut(&mut self) -> &mut R {
        self.queue.clear();
        self.is_eof = false;
        self.inner.as_mut().unwrap()
    }

    pub fn into_inner(mut self) -> R {
        self.inner.take().unwrap()
    }

    pub fn next_block(&mut self) -> io::Result<Option<Block>> {
        self.fill_queue()?;

        if let Some(buffered_rx) = self.queue.pop_front() {
            if let Ok(result) = buffered_rx.recv() {
                result.map(Some)
            } else {
                unreachable!();
            }
        } else {
            Ok(None)
        }
    }

    fn fill_queue(&mut self) -> io::Result<()> {
        use super::read_frame;

        let reader = self.inner.as_mut().unwrap();

        while self.queue.len() < self.queue.capacity() && !self.is_eof {
            match read_frame(reader)? {
                Some(buf) => {
                    let (buffered_tx, buffered_rx) = crossbeam_channel::bounded(1);

                    self.inflater_tx
                        .as_ref()
                        .unwrap()
                        .send((buf, buffered_tx))
                        .unwrap();

                    self.queue.push_back(buffered_rx);
                }
                None => self.is_eof = true,
            }
        }

        Ok(())
    }
}

fn spawn_inflaters(worker_count: usize, inflater_rx: InflaterRx) -> Vec<JoinHandle<()>> {
    use super::parse_frame;

    let mut handles = Vec::with_capacity(worker_count);

    for _ in 0..worker_count {
        let inflater_rx = inflater_rx.clone();

        handles.push(thread::spawn(move || {
            while let Ok((src, buffered_tx)) = inflater_rx.recv() {
                let result = parse_frame(&src);

                if buffered_tx.send(result).is_err() {
                    continue;
                }
            }
        }))
    }

    handles
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_get_mut() -> Result<(), Box<dyn std::error::Error>> {
        use crate::writer::BGZF_EOF;

        let worker_count = NonZeroUsize::try_from(2)?;
        let mut reader = Reader::with_worker_count(worker_count, BGZF_EOF);

        reader.fill_queue()?;

        assert!(!reader.queue.is_empty());
        assert!(reader.is_eof);

        let _ = reader.get_mut();

        assert!(reader.queue.is_empty());
        assert!(!reader.is_eof);

        Ok(())
    }
}
