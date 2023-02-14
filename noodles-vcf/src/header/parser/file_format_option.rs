use crate::header::parser::FileFormat;

/// A VCF header parser file format option.
#[derive(Debug, Default, Eq, PartialEq)]
pub enum FileFormatOption {
    /// Use the file format defined in the header.
    #[default]
    Auto,
    /// Override the file format with the given version.
    FileFormat(FileFormat),
}
