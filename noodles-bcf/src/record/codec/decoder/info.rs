mod field;

use std::io;

use noodles_vcf as vcf;

pub(crate) use self::field::read_field;
use crate::header::string_maps::StringStringMap;

pub fn read_info(
    src: &mut &[u8],
    infos: &vcf::header::Infos,
    string_string_map: &StringStringMap,
    len: usize,
) -> io::Result<vcf::record::Info> {
    let mut info = vcf::record::Info::default();

    for _ in 0..len {
        let (key, value) = read_field(src, infos, string_string_map)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        if info.insert(key.clone(), value).is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                vcf::record::info::TryFromFieldsError::DuplicateKey(key),
            ));
        }
    }

    Ok(info)
}
