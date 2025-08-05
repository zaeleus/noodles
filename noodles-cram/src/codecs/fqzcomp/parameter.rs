mod flags;

pub use self::flags::Flags;

use std::io::{self, Read};

use super::parameters::read_array;
use crate::io::reader::num::{read_u8, read_u16_le};

pub struct Parameter {
    pub context: u16,
    pub flags: Flags,
    pub max_sym: u8,
    pub q_bits: u8,
    pub q_shift: u8,
    pub q_loc: u8,
    pub s_loc: u8,
    pub p_loc: u8,
    pub d_loc: u8,
    pub q_map: Option<Vec<u8>>,
    pub q_tab: Vec<u8>,
    pub p_tab: Option<Vec<u8>>,
    pub d_tab: Option<Vec<u8>>,
}

pub fn fqz_decode_single_param<R>(reader: &mut R) -> io::Result<Parameter>
where
    R: Read,
{
    let context = read_u16_le(reader)?;
    let flags = read_u8(reader).map(Flags::from)?;
    let max_sym = read_u8(reader)?;

    let n = read_u8(reader)?;
    let (q_bits, q_shift) = (n >> 4, n & 0x0f);

    let n = read_u8(reader)?;
    let (q_loc, s_loc) = (n >> 4, n & 0x0f);

    let n = read_u8(reader)?;
    let (p_loc, d_loc) = (n >> 4, n & 0x0f);

    let q_map = if flags.contains(Flags::HAVE_QMAP) {
        let mut map = vec![0; usize::from(max_sym)];
        reader.read_exact(&mut map)?;
        Some(map)
    } else {
        None
    };

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
        read_array(reader, 1024).map(Some)?
    } else {
        None
    };

    let d_tab = if flags.contains(Flags::HAVE_DTAB) {
        read_array(reader, 256).map(Some)?
    } else {
        None
    };

    Ok(Parameter {
        context,
        flags,
        max_sym,
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
