pub mod op;

use std::str::FromStr;

pub use self::op::Op;

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

impl FromStr for Cigar {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut ops = Vec::new();

        let matches = s.match_indices(|c: char| !c.is_digit(10));
        let mut start = 0;

        for (end, raw_kind) in matches {
            let op = s[start..=end].parse().map_err(|_| ())?;
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
    fn test_from_str() {
        let actual = "1M13N144S".parse::<Cigar>().unwrap();

        let expected_ops = vec![
            Op::new(Kind::Match, 1),
            Op::new(Kind::Skip, 13),
            Op::new(Kind::SoftClip, 144),
        ];

        assert_eq!(actual.ops(), &expected_ops[..]);
    }
}
