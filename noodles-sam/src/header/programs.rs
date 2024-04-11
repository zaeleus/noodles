use bstr::BString;
use indexmap::IndexMap;

use super::record::value::{map::Program, Map};

type Inner = IndexMap<BString, Map<Program>>;

/// SAM header programs.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Programs(Inner);

impl AsRef<Inner> for Programs {
    fn as_ref(&self) -> &Inner {
        &self.0
    }
}

impl AsMut<Inner> for Programs {
    fn as_mut(&mut self) -> &mut Inner {
        &mut self.0
    }
}
