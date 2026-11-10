mod decode;
mod encode;

pub use self::{decode::decode, encode::encode};

const TYPE_COUNT: usize = 13;

#[allow(clippy::enum_variant_names)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Type {
    Type,
    String,
    Char,
    Digits0,
    DZLen,
    Dup,
    Diff,
    Digits,
    Delta,
    Delta0,
    Match,
    Nop,
    End,
}

impl From<Type> for u8 {
    fn from(ty: Type) -> Self {
        match ty {
            Type::Type => 0,
            Type::String => 1,
            Type::Char => 2,
            Type::Digits0 => 3,
            Type::DZLen => 4,
            Type::Dup => 5,
            Type::Diff => 6,
            Type::Digits => 7,
            Type::Delta => 8,
            Type::Delta0 => 9,
            Type::Match => 10,
            Type::Nop => 11,
            Type::End => 12,
        }
    }
}

impl From<Type> for usize {
    fn from(ty: Type) -> Self {
        Self::from(u8::from(ty))
    }
}

#[cfg(test)]
mod tests {
    use std::io;

    use super::*;

    #[test]
    fn test_self() -> io::Result<()> {
        let src = b"\
I17_08765:2:123:61541:01763#9\0\
I17_08765:2:123:1636:08611#9\0\
I17_08765:2:124:45613:16161#9\0\
";

        let input = encode(src)?;
        let output = decode(&input)?;

        assert_eq!(output, src);

        Ok(())
    }
}
