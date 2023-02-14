use crate::header::parser::FileFormat;

#[derive(Debug, Default, Eq, PartialEq)]
pub enum FileFormatOption {
    #[default]
    Auto,
    FileFormat(FileFormat),
}
