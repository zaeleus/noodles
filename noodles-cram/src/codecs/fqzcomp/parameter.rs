mod flags;

pub use self::flags::Flags;

use std::io::{self, Read};

use byteorder::{LittleEndian, ReadBytesExt};

use super::read_array;

pub struct Parameter {
    pub context: u16,
    pub flags: Flags,
    pub max_sym: u8,
    pub first_len: usize,
    pub last_len: usize,
    pub q_bits: u8,
    pub q_shift: u8,
    pub q_loc: u8,
    pub s_loc: u8,
    pub p_loc: u8,
    pub d_loc: u8,
    pub q_map: Vec<u8>,
    pub q_tab: Vec<u8>,
    pub p_tab: Option<Vec<u8>>,
    pub d_tab: Option<Vec<u8>>,
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

    if flags.contains(Flags::HAVE_QMAP) {
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
        first_len,
        last_len: 0,
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
