use std::io::{self, Write};

use crate::record::Chromosome;

pub(super) fn write_chromosome<W>(writer: &mut W, chromosome: &Chromosome) -> io::Result<()>
where
    W: Write,
{
    match chromosome {
        Chromosome::Name(name) => writer.write_all(name.as_bytes())?,
        Chromosome::Symbol(symbol) => {
            writer.write_all(b"<")?;
            writer.write_all(symbol.as_bytes())?;
            writer.write_all(b">")?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_chromosome() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, chromosome: &Chromosome, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_chromosome(buf, chromosome)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let chromosome = "sq0".parse()?;
        t(&mut buf, &chromosome, b"sq0")?;

        let chromosome = "<sq0>".parse()?;
        t(&mut buf, &chromosome, b"<sq0>")?;

        Ok(())
    }
}
