use super::{FileFormatOption, Parser};

/// A VCF header parser builder.
#[derive(Default)]
pub struct Builder {
    file_format_option: FileFormatOption,
}

impl Builder {
    /// Sets the file format option.
    pub fn set_file_format_option(mut self, file_format_option: FileFormatOption) -> Self {
        self.file_format_option = file_format_option;
        self
    }

    /// Builds a VCF header parser.
    pub fn build(self) -> Parser {
        Parser {
            file_format_option: self.file_format_option,
            ..Default::default()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();
        assert_eq!(builder.file_format_option, FileFormatOption::default());
    }
}
