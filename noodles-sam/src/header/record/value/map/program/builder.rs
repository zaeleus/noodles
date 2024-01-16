use super::Program;
use crate::header::record::value::map::{self, builder::BuildError};

/// A SAM header program builder.
#[derive(Debug, Default)]
pub struct Builder;

impl map::builder::Inner<Program> for Builder {
    fn build(self) -> Result<Program, BuildError> {
        Ok(Program)
    }
}
