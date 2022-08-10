use std::{
    collections::VecDeque,
    io::{self, Read},
    thread::{self, JoinHandle},
};

use bytes::Buf;
use crossbeam_channel::{Receiver, Sender};
use flate2::Crc;

use crate::{gz, Block, BGZF_HEADER_SIZE};

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
    pub(crate) fn with_worker_count(worker_count: usize, inner: R) -> Self {
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

    pub fn new(inner: R) -> Self {
        let worker_count = num_cpus::get();
        Self::with_worker_count(worker_count, inner)
    }

    pub fn get_ref(&self) -> &R {
        self.inner.as_ref().unwrap()
    }

    pub fn get_mut(&mut self) -> &mut R {
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

fn read_frame<R>(reader: &mut R) -> io::Result<Option<Vec<u8>>>
where
    R: Read,
{
    const MIN_FRAME_SIZE: usize = BGZF_HEADER_SIZE + gz::TRAILER_SIZE;
    const BSIZE_POSITION: usize = 16;

    let mut buf = vec![0; BGZF_HEADER_SIZE];

    match reader.read_exact(&mut buf) {
        Ok(()) => {}
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    }

    let bsize = (&buf[BSIZE_POSITION..]).get_u16_le();
    let block_size = usize::from(bsize) + 1;

    if block_size < MIN_FRAME_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid frame size",
        ));
    }

    buf.resize(block_size, 0);
    reader.read_exact(&mut buf[BGZF_HEADER_SIZE..])?;

    Ok(Some(buf))
}

fn split_frame(buf: &[u8]) -> (&[u8], &[u8], &[u8]) {
    let header = &buf[..BGZF_HEADER_SIZE];

    let n = buf.len() - gz::TRAILER_SIZE;
    let cdata = &buf[BGZF_HEADER_SIZE..n];

    let trailer = &buf[n..];

    (header, cdata, trailer)
}

fn parse_header(src: &[u8]) -> io::Result<()> {
    if is_valid_header(src) {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BGZF header",
        ))
    }
}

fn is_valid_header<B>(mut src: B) -> bool
where
    B: Buf,
{
    use std::mem;

    const BGZF_CM: u8 = 0x08; // DEFLATE
    const BGZF_FLG: u8 = 0x04; // FEXTRA
    const BGZF_XLEN: u16 = 6;
    const BGZF_SI1: u8 = b'B';
    const BGZF_SI2: u8 = b'C';
    const BGZF_SLEN: u16 = 2;

    let id_1 = src.get_u8();
    let id_2 = src.get_u8();
    let cm = src.get_u8();
    let flg = src.get_u8();

    // 4 (MTIME) + 1 (XFL) + 1 (OS)
    src.advance(mem::size_of::<u32>() + mem::size_of::<u8>() + mem::size_of::<u8>());

    let xlen = src.get_u16_le();
    let subfield_id_1 = src.get_u8();
    let subfield_id_2 = src.get_u8();
    let subfield_len = src.get_u16_le();

    id_1 == gz::MAGIC_NUMBER[0]
        && id_2 == gz::MAGIC_NUMBER[1]
        && cm == BGZF_CM
        && flg == BGZF_FLG
        && xlen == BGZF_XLEN
        && subfield_id_1 == BGZF_SI1
        && subfield_id_2 == BGZF_SI2
        && subfield_len == BGZF_SLEN
}

fn parse_trailer<B>(mut src: B) -> io::Result<(u32, usize)>
where
    B: Buf,
{
    let crc32 = src.get_u32_le();

    let r#isize = usize::try_from(src.get_u32_le())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((crc32, r#isize))
}

pub(crate) fn parse_frame(src: &[u8]) -> io::Result<Block> {
    let (header, cdata, trailer) = split_frame(src);

    parse_header(header)?;
    let (crc32, r#isize) = parse_trailer(trailer)?;

    let mut block = Block::default();

    let block_size =
        u64::try_from(src.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    block.set_size(block_size);

    let data = block.data_mut();
    data.set_position(0);
    data.resize(r#isize);

    inflate(cdata, crc32, data.as_mut())?;

    Ok(block)
}

fn inflate(src: &[u8], crc32: u32, dst: &mut [u8]) -> io::Result<()> {
    use super::inflate_data;

    inflate_data(src, dst)?;

    let mut crc = Crc::new();
    crc.update(dst);

    if crc.sum() == crc32 {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "block data checksum mismatch",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::writer::BGZF_EOF;

    #[test]
    fn test_parse_header() -> io::Result<()> {
        parse_header(BGZF_EOF)?;
        Ok(())
    }

    #[test]
    fn test_is_valid_header() {
        let mut src = [
            0x1f, 0x8b, // ID1, ID2
            0x08, // CM = DEFLATE
            0x04, // FLG = FEXTRA
            0x00, 0x00, 0x00, 0x00, // MTIME = 0
            0x00, // XFL = 0
            0xff, // OS = 255 (unknown)
            0x06, 0x00, // XLEN = 6
            b'B', b'C', // SI1, SI2
            0x02, 0x00, // SLEN = 2
            0x1b, 0x00, // BSIZE = 27
        ];

        let mut reader = &src[..];
        assert!(is_valid_header(&mut reader));

        src[0] = 0x00;
        let mut reader = &src[..];
        assert!(!is_valid_header(&mut reader));
    }

    #[test]
    fn test_parse_trailer() -> io::Result<()> {
        let (_, mut src) = BGZF_EOF.split_at(BGZF_EOF.len() - gz::TRAILER_SIZE);

        let (crc32, r#isize) = parse_trailer(&mut src)?;
        assert_eq!(crc32, 0);
        assert_eq!(r#isize, 0);

        Ok(())
    }

    #[test]
    fn test_read_frame() -> Result<(), Box<dyn std::error::Error>> {
        let mut src = BGZF_EOF;
        let buf = read_frame(&mut src)?.ok_or("invalid frame")?;
        assert_eq!(buf, BGZF_EOF);
        Ok(())
    }

    #[test]
    fn test_read_frame_with_invalid_block_size() {
        let data = {
            let mut eof = BGZF_EOF.to_vec();
            // BSIZE = 0
            eof[16] = 0x00;
            eof[17] = 0x00;
            eof
        };

        let mut reader = &data[..];
        assert!(read_frame(&mut reader).is_err());
    }
}
