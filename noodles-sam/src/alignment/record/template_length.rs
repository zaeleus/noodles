use std::io;

/// An alignment record template length.
pub trait TemplateLength {
    /// Converts a mapping quality to a `i32`.
    fn try_to_i32(&self) -> io::Result<i32>;
}

impl TryFrom<&dyn TemplateLength> for i32 {
    type Error = io::Error;

    fn try_from(raw_template_length: &dyn TemplateLength) -> Result<Self, Self::Error> {
        raw_template_length.try_to_i32()
    }
}
