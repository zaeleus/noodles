use std::io;

/// Raw SAM record template length.
#[derive(Debug, Eq, PartialEq)]
pub struct TemplateLength<'a>(&'a [u8]);

impl<'a> TemplateLength<'a> {
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self(buf)
    }
}

impl<'a> AsRef<[u8]> for TemplateLength<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> TryFrom<TemplateLength<'a>> for i32 {
    type Error = io::Error;

    fn try_from(raw_template_length: TemplateLength<'a>) -> Result<Self, Self::Error> {
        lexical_core::parse(raw_template_length.as_ref())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}
