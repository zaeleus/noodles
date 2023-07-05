#![allow(dead_code)]

mod group;

use std::ops::Deref;

use self::group::Group;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BaseModifications(Vec<Group>);

impl Deref for BaseModifications {
    type Target = [Group];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
