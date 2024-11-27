use std::io::{self, Write};

use crate::directive_buf::value::GenomeBuild;

pub(crate) fn write_genome_build<W>(writer: &mut W, genome_build: &GenomeBuild) -> io::Result<()>
where
    W: Write,
{
    write!(writer, "{} {}", genome_build.source(), genome_build.name())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_genome_build() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();
        let genome_build = GenomeBuild::new(String::from("NDLS"), String::from("r1"));
        write_genome_build(&mut buf, &genome_build)?;
        assert_eq!(buf, b"NDLS r1");
        Ok(())
    }
}
