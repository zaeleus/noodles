mod flags;

use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use self::flags::Flags;
use super::read_array;

pub struct Parameter {
    context: u16,
    flags: Flags,
    pub max_sym: u8,
    first_len: usize,
    q_bits: u8,
    q_shift: u8,
    q_loc: u8,
    s_loc: u8,
    p_loc: u8,
    d_loc: u8,
    q_map: Vec<u8>,
    q_tab: Vec<u8>,
    p_tab: Vec<u8>,
    d_tab: Vec<u8>,
}

pub fn fqz_decode_single_param<R>(reader: &mut R) -> io::Result<Parameter>
where
    R: Read,
{
    let context = reader.read_u16::<LittleEndian>()?;
    let flags = reader.read_u8().map(Flags::from)?;
    let max_sym = reader.read_u8()?;
    let first_len = 1;

    let x = reader.read_u8()?;
    let q_bits = x / 16;
    let q_shift = x % 16;

    let x = reader.read_u8()?;
    let q_loc = x / 16;
    let s_loc = x % 16;

    let x = reader.read_u8()?;
    let p_loc = x / 16;
    let d_loc = x % 16;

    let mut q_map = vec![0; usize::from(max_sym)];

    if flags.contains(Flags::DO_DEDUP) {
        reader.read_exact(&mut q_map)?;
    }

    let q_tab = if flags.contains(Flags::HAVE_QTAB) {
        read_array(reader, 256)?
    } else {
        let mut tab = Vec::with_capacity(256);

        for i in 0..=u8::MAX {
            tab.push(i);
        }

        tab
    };

    let p_tab = if flags.contains(Flags::HAVE_PTAB) {
        read_array(reader, 1024)?
    } else {
        vec![0; 1024]
    };

    let d_tab = if flags.contains(Flags::HAVE_DTAB) {
        read_array(reader, 256)?
    } else {
        vec![0; 256]
    };

    Ok(Parameter {
        context,
        flags,
        max_sym,
        first_len,
        q_bits,
        q_shift,
        q_loc,
        s_loc,
        p_loc,
        d_loc,
        q_map,
        q_tab,
        p_tab,
        d_tab,
    })
}
