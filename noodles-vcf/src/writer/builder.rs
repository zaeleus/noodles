use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use noodles_bgzf as bgzf;

use super::Writer;

/// A BAM writer builder.
#[derive(Debug, Default)]
pub struct Builder;

impl Builder {
    /// Builds a VCF writer from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_vcf as vcf;
    /// let writer = vcf::writer::Builder.build_from_path("out.vcf")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, dst: P) -> io::Result<Writer<Box<dyn Write>>>
    where
        P: AsRef<Path>,
    {
        let dst = dst.as_ref();

        let file = File::create(dst)?;

        let writer: Box<dyn Write> = match dst.extension().and_then(|ext| ext.to_str()) {
            Some("gz" | "bgz") => Box::new(bgzf::Writer::new(file)),
            _ => Box::new(BufWriter::new(file)),
        };

        Ok(Writer::new(writer))
    }
}
