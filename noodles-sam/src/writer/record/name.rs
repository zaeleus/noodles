use std::io::{self, Write};

use crate::record::Name;

pub(super) fn write_name<W>(writer: &mut W, name: Option<&Name>) -> io::Result<()>
where
    W: Write,
{
    const MISSING: &[u8] = b"*";

    let qname = name.map(|s| s.as_ref()).unwrap_or(MISSING);
    writer.write_all(qname)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_name() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        buf.clear();
        write_name(&mut buf, None)?;
        assert_eq!(buf, b"*");

        buf.clear();
        let name = Name::try_new(b"r0".to_vec())?;
        write_name(&mut buf, Some(&name))?;
        assert_eq!(buf, b"r0");

        Ok(())
    }
}
