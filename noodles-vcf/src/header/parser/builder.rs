use super::{FileFormatOption, Parser};

#[derive(Default)]
pub struct Builder {
    file_format_option: FileFormatOption,
}

impl Builder {
    pub fn set_file_format_option(mut self, file_format_option: FileFormatOption) -> Self {
        self.file_format_option = file_format_option;
        self
    }

    pub fn build(self) -> Parser {
        Parser {
            file_format_option: self.file_format_option,
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
