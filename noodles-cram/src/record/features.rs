use std::ops::{Deref, DerefMut};

use super::Feature;

/// CRAM record features.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Features(Vec<Feature>);

impl Deref for Features {
    type Target = Vec<Feature>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Features {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl From<Vec<Feature>> for Features {
    fn from(features: Vec<Feature>) -> Self {
        Self(features)
    }
}
