mod file_id;
mod format_version;
pub(crate) mod magic_number;

use std::io::{self, Read};

use self::{
    file_id::read_file_id, format_version::read_format_version, magic_number::read_magic_number,
};
use crate::FileDefinition;

pub(super) fn read_file_definition<R>(reader: &mut R) -> io::Result<FileDefinition>
where
    R: Read,
{
    read_magic_number(reader).and_then(magic_number::validate)?;

    let version = read_format_version(reader)?;
    let file_id = read_file_id(reader)?;

    Ok(FileDefinition::new(version, file_id))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::file_definition::Version;

    #[test]
    fn test_read_file_definition() -> Result<(), Box<dyn std::error::Error>> {
        let src = [
            0x43, 0x52, 0x41, 0x4d, // magic number = b"CRAM"
            0x03, 0x00, // format version = (3, 0)
            0x00, 0x68, 0xac, 0xf3, 0x06, 0x4d, 0xaa, 0x1e, 0x29, 0xa4, 0xa0, 0x8c, 0x56, 0xee,
            0x91, 0x9b, 0x91, 0x04, 0x21, 0x1f, // file ID
        ];

        let mut reader = &src[..];
        let actual = read_file_definition(&mut reader)?;

        let file_id = <[u8; 20]>::try_from(&src[6..])?;
        let expected = FileDefinition::new(Version::new(3, 0), file_id);

        assert_eq!(actual, expected);

        Ok(())
    }
}
