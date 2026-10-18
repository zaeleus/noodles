use std::{
    error, fmt,
    io::{self, Write},
};

use crate::header::FileFormat;

const VCF_4_3: FileFormat = FileFormat::new(4, 3);
const VCF_4_4: FileFormat = FileFormat::new(4, 4);

/// An error returned when a reference sequence name fails to write.
#[derive(Debug)]
pub enum WriteError {
    /// I/O error.
    Io(io::Error),
    /// The reference sequence name is invalid.
    InvalidReferenceSequenceName(String),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidReferenceSequenceName(_) => None,
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidReferenceSequenceName(s) => {
                write!(f, "invalid reference sequence name: {s}")
            }
        }
    }
}

pub(super) fn write_reference_sequence_name<W>(
    writer: &mut W,
    file_format: FileFormat,
    reference_sequence_name: &str,
) -> Result<(), WriteError>
where
    W: Write,
{
    let name = strip_symbol_delimiters(reference_sequence_name).unwrap_or(reference_sequence_name);

    if is_valid(file_format, name) {
        writer
            .write_all(reference_sequence_name.as_bytes())
            .map_err(WriteError::Io)
    } else {
        Err(WriteError::InvalidReferenceSequenceName(
            reference_sequence_name.into(),
        ))
    }
}

fn strip_symbol_delimiters(s: &str) -> Option<&str> {
    const PREFIX: char = '<';
    const SUFFIX: char = '>';

    s.strip_prefix(PREFIX).and_then(|t| t.strip_suffix(SUFFIX))
}

fn is_valid(file_format: FileFormat, s: &str) -> bool {
    if file_format < VCF_4_3 {
        !s.is_empty() && !s.chars().any(|c| c.is_whitespace())
    } else {
        let mut chars = s.chars();

        let is_valid_first_char = chars
            .next()
            .map(|c| c != '*' && c != '=' && is_valid_char(c))
            .unwrap_or_default();

        is_valid_first_char
            && chars.all(is_valid_char)
            && !is_reserved_symbolic_allele_name(file_format, s)
    }
}

fn is_reserved_symbolic_allele_name(file_format: FileFormat, s: &str) -> bool {
    let mut is_reserved = matches!(
        s,
        "DEL" | "INS" | "DUP" | "INV" | "CNV" | "DUP:TANDEM" | "DEL:ME" | "INS:ME" | "NON_REF"
    );

    is_reserved = if file_format < VCF_4_4 {
        is_reserved || s == "BND"
    } else {
        is_reserved || s == "CNV:TR"
    };

    is_reserved
}

// § 1.4.7 "Contig field format" (2025-11-05): "`[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*`".
fn is_valid_char(c: char) -> bool {
    ('!'..='~').contains(&c)
        && !matches!(
            c,
            '\\' | ',' | '"' | '`' | '\'' | '(' | ')' | '[' | ']' | '{' | '}' | '<' | '>',
        )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_reference_sequence_name() -> Result<(), WriteError> {
        fn t(
            buf: &mut Vec<u8>,
            file_format: FileFormat,
            reference_sequence_name: &str,
            expected: &[u8],
        ) -> Result<(), WriteError> {
            buf.clear();
            write_reference_sequence_name(buf, file_format, reference_sequence_name)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, VCF_4_3, "sq0", b"sq0")?;
        t(&mut buf, VCF_4_3, "<sq0>", b"<sq0>")?;

        buf.clear();
        assert!(matches!(
            write_reference_sequence_name(&mut buf, VCF_4_3, "sq 0"),
            Err(WriteError::InvalidReferenceSequenceName(s)) if s == "sq 0"
        ));

        buf.clear();
        assert!(matches!(
            write_reference_sequence_name(&mut buf, VCF_4_3, "<sq 0>"),
            Err(WriteError::InvalidReferenceSequenceName(s)) if s == "<sq 0>"
        ));

        Ok(())
    }

    #[test]
    fn test_strip_symbol_delimiters() {
        assert!(strip_symbol_delimiters("sq0").is_none());
        assert_eq!(strip_symbol_delimiters("<sq0>"), Some("sq0"));
    }

    #[test]
    fn test_is_valid() {
        const VCF_4_2: FileFormat = FileFormat::new(4, 2);

        assert!(is_valid(VCF_4_2, "sq0"));
        assert!(!is_valid(VCF_4_2, ""));
        assert!(!is_valid(VCF_4_2, "sq 0"));
        assert!(is_valid(VCF_4_2, "sq[0]"));
        assert!(is_valid(VCF_4_2, ">sq0"));
        assert!(is_valid(VCF_4_2, "*sq0"));
        assert!(is_valid(VCF_4_2, "=sq0"));
        assert!(!is_valid(VCF_4_2, "\t"));
        assert!(is_valid(VCF_4_2, "DEL"));

        assert!(is_valid(VCF_4_3, "sq0"));
        assert!(!is_valid(VCF_4_3, ""));
        assert!(!is_valid(VCF_4_3, "sq 0"));
        assert!(!is_valid(VCF_4_3, "sq[0]"));
        assert!(!is_valid(VCF_4_3, ">sq0"));
        assert!(!is_valid(VCF_4_3, "*sq0"));
        assert!(!is_valid(VCF_4_3, "=sq0"));
        assert!(!is_valid(VCF_4_3, "\t"));
        assert!(!is_valid(VCF_4_3, "DEL"));
    }

    #[test]
    fn test_is_reserved_symbolic_allele_name() {
        assert!(is_reserved_symbolic_allele_name(VCF_4_3, "DEL"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_3, "INS"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_3, "DUP"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_3, "INV"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_3, "CNV"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_3, "BND"));
        assert!(!is_reserved_symbolic_allele_name(VCF_4_3, "CNV:TR"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_3, "DUP:TANDEM"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_3, "DEL:ME"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_3, "INS:ME"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_3, "NON_REF"));

        assert!(is_reserved_symbolic_allele_name(VCF_4_4, "DEL"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_4, "INS"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_4, "DUP"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_4, "INV"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_4, "CNV"));
        assert!(!is_reserved_symbolic_allele_name(VCF_4_4, "BND"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_4, "CNV:TR"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_4, "DUP:TANDEM"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_4, "DEL:ME"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_4, "INS:ME"));
        assert!(is_reserved_symbolic_allele_name(VCF_4_4, "NON_REF"));
    }
}
