use std::io;

use noodles_sam as sam;

/// A raw BAM record template length.
#[derive(Debug, Eq, PartialEq)]
pub struct TemplateLength(i32);

impl TemplateLength {
    pub(super) fn new(n: i32) -> Self {
        Self(n)
    }
}

impl sam::alignment::record::field::TemplateLength for TemplateLength {
    fn try_to_i32(&self) -> io::Result<i32> {
        Ok(self.0)
    }
}

impl From<TemplateLength> for i32 {
    fn from(template_length: TemplateLength) -> Self {
        template_length.0
    }
}
