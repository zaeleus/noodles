use std::io;

use byteorder::{ByteOrder, LittleEndian};
use bytes::BytesMut;
use tokio_util::codec::Decoder;

use crate::BGZF_HEADER_SIZE;

pub struct BlockCodec;

impl Decoder for BlockCodec {
    type Item = BytesMut;
    type Error = io::Error;

    fn decode(&mut self, src: &mut BytesMut) -> Result<Option<Self::Item>, Self::Error> {
        if src.len() < BGZF_HEADER_SIZE {
            src.reserve(BGZF_HEADER_SIZE);
            return Ok(None);
        }

        let block_size = {
            let header = &src[..BGZF_HEADER_SIZE];
            let bsize = LittleEndian::read_u16(&header[16..]);
            usize::from(bsize) + 1
        };

        if src.len() < block_size {
            src.reserve(block_size);
            return Ok(None);
        }

        Ok(Some(src.split_to(block_size)))
    }
}
