/// Raw SAM record flags.
#[derive(Debug, Eq, PartialEq)]
pub struct Flags<'a>(&'a [u8]);

impl<'a> Flags<'a> {
    pub(super) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }
}

impl<'a> AsRef<[u8]> for Flags<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> TryFrom<Flags<'a>> for u16 {
    type Error = lexical_core::Error;

    fn try_from(flags: Flags<'a>) -> Result<Self, Self::Error> {
        lexical_core::parse(flags.as_ref())
    }
}

impl<'a> TryFrom<Flags<'a>> for crate::record::Flags {
    type Error = lexical_core::Error;

    fn try_from(flags: Flags<'a>) -> Result<Self, Self::Error> {
        u16::try_from(flags).map(Self::from)
    }
}
