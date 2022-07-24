use super::Meta;
use crate::header::record::value::map::{self, builder::BuildError};

#[derive(Default)]
pub struct Builder {
    values: Option<Vec<String>>,
}

impl map::builder::Inner<Meta> for Builder {
    fn build(self) -> Result<Meta, BuildError> {
        let values = self.values.ok_or(BuildError::MissingField("Values"))?;
        Ok(Meta { values })
    }
}

impl map::Builder<Meta> {
    /// Sets the values.
    pub fn set_values(mut self, values: Vec<String>) -> Self {
        self.inner.values = Some(values);
        self
    }
}
