use super::AlternativeAllele;
use crate::header::record::value::map::{self, builder::BuildError};

#[derive(Default)]
pub struct Builder {
    description: Option<String>,
}

impl map::builder::Inner<AlternativeAllele> for Builder {
    fn build(self) -> Result<AlternativeAllele, BuildError> {
        let description = self.description.ok_or(BuildError::MissingField("Values"))?;
        Ok(AlternativeAllele { description })
    }
}

impl map::builder::Described<AlternativeAllele> for Builder {
    fn set_description<D>(mut self, description: D) -> Self
    where
        D: Into<String>,
    {
        self.description = Some(description.into());
        self
    }
}
