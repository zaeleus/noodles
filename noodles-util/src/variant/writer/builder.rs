use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf as vcf;

use super::Writer;
use crate::variant::{CompressionMethod, Format};

/// A variant writer builder.
#[derive(Default)]
pub struct Builder {
    compression_method: Option<Option<CompressionMethod>>,
    format: Option<Format>,
}

impl Builder {
    /// Sets the compression method of the output.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{self, CompressionMethod};
    /// let builder = variant::writer::Builder::default()
    ///     .set_compression_method(Some(CompressionMethod::Bgzf));
    /// ```
    pub fn set_compression_method(mut self, compression_method: Option<CompressionMethod>) -> Self {
        self.compression_method = Some(compression_method);
        self
    }

    /// Sets the format of the output.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{self, Format};
    /// let builder = variant::writer::Builder::default().set_format(Format::Vcf);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Builds a variant writer from a path.
    ///
    /// If the format or compression method is not set, it is detected from the path extension.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_util::variant::{self, Format};
    /// let writer = variant::writer::Builder::default().build_from_path("out.vcf.gz")?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_path<P>(mut self, src: P) -> io::Result<Writer>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        if self.compression_method.is_none() {
            self.compression_method = Some(detect_compression_method_from_path_extension(src));
        }

        if self.format.is_none() {
            self.format = detect_format_from_path_extension(src);
        }

        let file = File::create(src).map(BufWriter::new)?;
        Ok(self.build_from_writer(file))
    }

    /// Builds a variant writer from a writer.
    ///
    /// If the format is not set, a default format is used. If the compression method is not set, a
    /// default one is determined by the format.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::variant;
    /// let writer = variant::writer::Builder::default().build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> Writer
    where
        W: Write + 'static,
    {
        let format = self.format.unwrap_or(Format::Vcf);

        let compression_method = match self.compression_method {
            Some(compression_method) => compression_method,
            None => match format {
                Format::Vcf => None,
                Format::Bcf => Some(CompressionMethod::Bgzf),
            },
        };

        let inner: Box<dyn vcf::VariantWriter> = match (format, compression_method) {
            (Format::Vcf, None) => Box::new(vcf::Writer::new(writer)),
            (Format::Vcf, Some(CompressionMethod::Bgzf)) => {
                Box::new(vcf::Writer::new(bgzf::Writer::new(writer)))
            }
            (Format::Bcf, None) => Box::new(bcf::Writer::from(writer)),
            (Format::Bcf, Some(CompressionMethod::Bgzf)) => Box::new(bcf::Writer::new(writer)),
        };

        Writer { inner }
    }
}

fn detect_format_from_path_extension<P>(path: P) -> Option<Format>
where
    P: AsRef<Path>,
{
    let path = path.as_ref();

    match path.extension().and_then(|ext| ext.to_str()) {
        Some("vcf") => Some(Format::Vcf),
        Some("bcf") => Some(Format::Bcf),
        Some("gz" | "bgz") => {
            let file_stem = path.file_stem().and_then(|stem| stem.to_str())?;
            file_stem.ends_with("vcf").then_some(Format::Vcf)
        }
        _ => None,
    }
}

fn detect_compression_method_from_path_extension<P>(path: P) -> Option<CompressionMethod>
where
    P: AsRef<Path>,
{
    match path.as_ref().extension().and_then(|ext| ext.to_str()) {
        Some("bcf" | "gz" | "bgz") => Some(CompressionMethod::Bgzf),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_compression_method_from_path_extension() {
        assert_eq!(
            detect_compression_method_from_path_extension("out.vcf"),
            None
        );
        assert_eq!(
            detect_compression_method_from_path_extension("out.vcf.gz"),
            Some(CompressionMethod::Bgzf),
        );
        assert_eq!(
            detect_compression_method_from_path_extension("out.bcf"),
            Some(CompressionMethod::Bgzf)
        );
    }

    #[test]
    fn test_detect_format_from_path_extension() {
        assert_eq!(
            detect_format_from_path_extension("out.vcf"),
            Some(Format::Vcf)
        );
        assert_eq!(
            detect_format_from_path_extension("out.vcf.gz"),
            Some(Format::Vcf)
        );
        assert_eq!(
            detect_format_from_path_extension("out.bcf"),
            Some(Format::Bcf)
        );

        assert!(detect_format_from_path_extension("out.bcf.gz").is_none());
        assert!(detect_format_from_path_extension("out.fa").is_none());
    }
}
