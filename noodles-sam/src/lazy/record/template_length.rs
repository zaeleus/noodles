use std::io;

use crate::alignment::record::TemplateLength as _;

/// Raw SAM record template length.
#[derive(Debug, Eq, PartialEq)]
pub struct TemplateLength<'a>(&'a [u8]);

impl<'a> TemplateLength<'a> {
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self(buf)
    }
}

impl<'a> crate::alignment::record::TemplateLength for TemplateLength<'a> {
    fn try_to_i32(&self) -> io::Result<i32> {
        lexical_core::parse(self.as_ref())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
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
        raw_template_length.try_to_i32()
    }
}
