use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

use noodles_bam as bam;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::IndexedReader;
use crate::alignment::Format;

/// An indexed alignment reader builder.
#[derive(Default)]
pub struct Builder {
    format: Option<Format>,
    reference_sequence_repository: fasta::Repository,
}

impl Builder {
    /// Sets the format of the input.
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Sets the reference sequence repository.
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Builds an indexed alignment reader from a path.
    pub fn build_from_path<P>(self, src: P) -> io::Result<IndexedReader<File>>
    where
        P: AsRef<Path>,
    {
        use crate::alignment::reader::builder::detect_format;

        let mut reader = File::open(src.as_ref()).map(BufReader::new)?;

        let format = match self.format {
            Some(format) => format,
            None => detect_format(&mut reader)?,
        };

        match format {
            Format::Sam => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "source not bgzip-compressed",
            )),
            Format::SamGz => sam::indexed_reader::Builder::default()
                .build_from_path(src)
                .map(IndexedReader::Sam),
            Format::Bam => bam::indexed_reader::Builder::default()
                .build_from_path(src)
                .map(IndexedReader::Bam),
            Format::Cram => cram::indexed_reader::Builder::default()
                .set_reference_sequence_repository(self.reference_sequence_repository)
                .build_from_path(src)
                .map(IndexedReader::Cram),
        }
    }
}
