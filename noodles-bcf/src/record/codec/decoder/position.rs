use std::io;

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_core::Position;

pub fn read_pos(src: &mut &[u8]) -> io::Result<Option<Position>> {
    const TELOMERE_START: i32 = -1;

    match src.read_i32::<LittleEndian>()? {
        TELOMERE_START => Ok(None),
        n => usize::try_from(n)
            .map(|m| m + 1)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .map(Position::new),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_pos() -> io::Result<()> {
        let data = (-1i32).to_le_bytes();
        let mut src = &data[..];
        assert!(read_pos(&mut src)?.is_none());

        let data = 0i32.to_le_bytes();
        let mut src = &data[..];
        assert_eq!(read_pos(&mut src)?, Some(Position::MIN));

        let data = (-2i32).to_le_bytes();
        let mut src = &data[..];
        assert!(matches!(
            read_pos(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[cfg(not(target_pointer_width = "16"))]
    #[test]
    fn test_read_pos_with_max_position() -> Result<(), Box<dyn std::error::Error>> {
        let data = i32::MAX.to_le_bytes();
        let mut src = &data[..];
        let expected = Position::try_from(1 << 31)?;
        assert_eq!(read_pos(&mut src)?, Some(expected));
        Ok(())
    }
}
