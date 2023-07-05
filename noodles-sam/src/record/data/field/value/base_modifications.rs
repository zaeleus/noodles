#![allow(dead_code)]

pub mod group;
mod parser;

pub use self::group::Group;

use std::ops::Deref;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BaseModifications(Vec<Group>);

impl Deref for BaseModifications {
    type Target = [Group];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
