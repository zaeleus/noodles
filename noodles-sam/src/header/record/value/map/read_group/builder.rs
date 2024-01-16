use super::ReadGroup;
use crate::header::record::value::map::{self, builder::BuildError};

/// A SAM header reference read group.
#[derive(Debug, Default)]
pub struct Builder;

impl map::builder::Inner<ReadGroup> for Builder {
    fn build(self) -> Result<ReadGroup, BuildError> {
        Ok(ReadGroup)
    }
}
