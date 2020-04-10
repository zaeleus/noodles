pub mod op;

use std::{fmt, str::FromStr};

pub use self::op::Op;

use super::record;

#[derive(Debug, Default)]
pub struct Cigar {
    ops: Vec<Op>,
}

impl Cigar {
    pub fn new(ops: Vec<Op>) -> Self {
        Self { ops }
    }

    pub fn ops(&self) -> &[Op] {
        &self.ops
    }
}

impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for op in self.ops() {
            write!(f, "{}", op)?;
        }

        Ok(())
    }
}

impl FromStr for Cigar {
    type Err = op::ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(op::ParseError::Empty);
        } else if s == record::NULL_FIELD {
            return Ok(Cigar::default());
        }

        let mut ops = Vec::new();

        let matches = s.match_indices(|c: char| !c.is_digit(10));
        let mut start = 0;

        for (end, raw_kind) in matches {
            let op = s[start..=end].parse()?;
            ops.push(op);
            start = end + raw_kind.len();
        }

        Ok(Cigar::new(ops))
    }
}

#[cfg(test)]
mod tests {
    use super::{op::Kind, *};

    #[test]
    fn test_fmt() {
        let cigar = Cigar::new(vec![
            Op::new(Kind::Match, 1),
            Op::new(Kind::Skip, 13),
            Op::new(Kind::SoftClip, 144),
        ]);

        assert_eq!(format!("{}", cigar), "1M13N144S");
    }

    #[test]
    fn test_from_str() {
        let actual = "1M13N144S".parse::<Cigar>().unwrap();

        let expected_ops = vec![
            Op::new(Kind::Match, 1),
            Op::new(Kind::Skip, 13),
            Op::new(Kind::SoftClip, 144),
        ];

        assert_eq!(actual.ops(), &expected_ops[..]);

        assert!("".parse::<Cigar>().is_err());
    }

    #[test]
    fn test_from_str_with_null_cigar() {
        let cigar = "*".parse::<Cigar>().unwrap();
        assert!(cigar.ops().is_empty());
    }
}
