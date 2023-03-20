use std::io;

use crate::record::Chromosome;

pub(super) fn parse_chromosome(s: &str, chromosome: &mut Chromosome) -> io::Result<()> {
    // symbol
    if let Some(t) = s.strip_prefix('<') {
        if let Some(t) = t.strip_suffix('>') {
            *chromosome = Chromosome::Symbol(t.into());
            return Ok(());
        }
    }

    // name
    if is_valid_name(s) {
        *chromosome = Chromosome::Name(s.into());
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid chromosome",
        ))
    }
}

// § 1.4.7 "Contig field format"
fn is_valid_name_char(c: char) -> bool {
    ('!'..='~').contains(&c)
        && !matches!(
            c,
            '\\' | ',' | '"' | '`' | '\'' | '(' | ')' | '[' | ']' | '{' | '}' | '<' | '>',
        )
}

fn is_valid_name(s: &str) -> bool {
    if s.is_empty() {
        return false;
    }

    let mut chars = s.chars();

    if let Some(c) = chars.next() {
        if c == '*' || c == '=' || !is_valid_name_char(c) {
            return false;
        }
    }

    chars.all(is_valid_name_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_chromosome() -> io::Result<()> {
        let mut chromosome = Chromosome::Name(String::from("."));

        parse_chromosome("sq0", &mut chromosome)?;
        assert_eq!(chromosome, Chromosome::Name(String::from("sq0")));

        parse_chromosome("<sq0>", &mut chromosome)?;
        assert_eq!(chromosome, Chromosome::Symbol(String::from("sq0")));

        assert!(matches!(
            parse_chromosome("", &mut chromosome),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid_name() {
        assert!(is_valid_name("sq0"));

        assert!(!is_valid_name(""));
        assert!(!is_valid_name("sq 0"));
        assert!(!is_valid_name("sq[0]"));
        assert!(!is_valid_name(">sq0"));
        assert!(!is_valid_name("*sq0"));
        assert!(!is_valid_name("=sq0"));
    }
}
