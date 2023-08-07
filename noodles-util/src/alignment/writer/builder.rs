//! Alignment writer builder.

use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use cram::data_container::BlockContentEncoderMap;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::Writer;
use crate::alignment::{CompressionMethod, Format};

/// An alignment writer builder.
#[derive(Default)]
pub struct Builder {
    compression_method: Option<Option<CompressionMethod>>,
    format: Option<Format>,
    reference_sequence_repository: fasta::Repository,
    block_content_encoder_map: BlockContentEncoderMap,
}

impl Builder {
    /// Sets the compression method.
    pub fn set_compression_method(mut self, compression_method: Option<CompressionMethod>) -> Self {
        self.compression_method = Some(compression_method);
        self
    }

    /// Sets the format of the output.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::alignment::{self, Format};
    /// let builder = alignment::writer::Builder::default().set_format(Format::Sam);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Sets the reference sequence repository.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let repository = fasta::Repository::default();
    ///
    /// let builder = alignment::writer::Builder::default()
    ///     .set_reference_sequence_repository(repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Sets the block content-encoder map.
    ///
    /// This is only used when the output format is CRAM.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram::data_container::BlockContentEncoderMap;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let builder = alignment::writer::Builder::default()
    ///     .set_format(Format::Cram)
    ///     .set_block_content_encoder_map(BlockContentEncoderMap::default());
    /// ```
    pub fn set_block_content_encoder_map(
        mut self,
        block_content_encoder_map: BlockContentEncoderMap,
    ) -> Self {
        self.block_content_encoder_map = block_content_encoder_map;
        self
    }

    /// Builds an alignment writer from a path.
    ///
    /// If the format is not set, it is detected from the path extension.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let writer = alignment::writer::Builder::default()
    ///     .set_format(Format::Sam)
    ///     .build_from_path("out.sam")?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_path<P>(mut self, src: P) -> io::Result<Writer>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        if self.format.is_none() {
            if let Some((format, compression_method)) =
                detect_format_and_compression_method_from_path_extension(src)
            {
                self.format = Some(format);
                self.compression_method = Some(compression_method);
            }
        }

        File::create(src)
            .map(BufWriter::new)
            .and_then(|writer| self.build_from_writer(writer))
    }

    /// Builds an alignment writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let writer = alignment::writer::Builder::default()
    ///     .set_format(Format::Sam)
    ///     .build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> io::Result<Writer>
    where
        W: Write + 'static,
    {
        let format = self.format.unwrap_or(Format::Sam);

        let compression_method = match self.compression_method {
            Some(compression_method) => compression_method,
            None => match format {
                Format::Sam | Format::Cram => None,
                Format::Bam => Some(CompressionMethod::Bgzf),
            },
        };

        let inner: Box<dyn sam::AlignmentWriter> = match (format, compression_method) {
            (Format::Sam, None) => Box::new(sam::Writer::new(writer)),
            (Format::Sam, Some(CompressionMethod::Bgzf)) => {
                Box::new(sam::Writer::new(bgzf::Writer::new(writer)))
            }
            (Format::Bam, None) => Box::new(bam::Writer::from(writer)),
            (Format::Bam, Some(CompressionMethod::Bgzf)) => Box::new(bam::Writer::new(writer)),
            (Format::Cram, None) => Box::new(
                cram::writer::Builder::default()
                    .set_reference_sequence_repository(self.reference_sequence_repository)
                    .set_block_content_encoder_map(self.block_content_encoder_map)
                    .build_with_writer(writer),
            ),
            (Format::Cram, Some(CompressionMethod::Bgzf)) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "CRAM cannot be bgzip-compressed",
                ));
            }
        };

        Ok(Writer { inner })
    }
}

fn detect_format_and_compression_method_from_path_extension<P>(
    path: P,
) -> Option<(Format, Option<CompressionMethod>)>
where
    P: AsRef<Path>,
{
    let ext = path.as_ref().file_name().and_then(|ext| ext.to_str())?;

    if ext.ends_with("sam") {
        Some((Format::Sam, None))
    } else if ext.ends_with("sam.gz") || ext.ends_with("sam.bgz") {
        Some((Format::Sam, Some(CompressionMethod::Bgzf)))
    } else if ext.ends_with("bam") {
        Some((Format::Bam, Some(CompressionMethod::Bgzf)))
    } else if ext.ends_with("cram") {
        Some((Format::Cram, None))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_format_and_compression_method_from_path_extension() {
        assert_eq!(
            detect_format_and_compression_method_from_path_extension("out.sam"),
            Some((Format::Sam, None))
        );
        assert_eq!(
            detect_format_and_compression_method_from_path_extension("out.sam.gz"),
            Some((Format::Sam, Some(CompressionMethod::Bgzf)))
        );
        assert_eq!(
            detect_format_and_compression_method_from_path_extension("out.bam"),
            Some((Format::Bam, Some(CompressionMethod::Bgzf)))
        );
        assert_eq!(
            detect_format_and_compression_method_from_path_extension("out.cram"),
            Some((Format::Cram, None))
        );

        assert!(detect_format_and_compression_method_from_path_extension("out.fa").is_none());
    }
}
