/// A raw BAM record template length.
#[derive(Debug, Eq, PartialEq)]
pub struct TemplateLength(i32);

impl TemplateLength {
    pub(super) fn new(n: i32) -> Self {
        Self(n)
    }
}

impl From<TemplateLength> for i32 {
    fn from(template_length: TemplateLength) -> Self {
        template_length.0
    }
}
