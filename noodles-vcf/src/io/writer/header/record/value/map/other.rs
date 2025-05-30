use std::io::{self, Write};

use super::write_other_fields;
use crate::header::record::value::{Map, map::Other};

pub(crate) fn write_other<W>(writer: &mut W, other: &Map<Other>) -> io::Result<()>
where
    W: Write,
{
    write_other_fields(writer, other.other_fields())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_other() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        buf.clear();
        let map = Map::<Other>::builder()
            .insert("noodles".parse()?, "vcf")
            .build()?;
        write_other(&mut buf, &map)?;
        assert_eq!(buf, br#",noodles="vcf""#);

        Ok(())
    }
}
