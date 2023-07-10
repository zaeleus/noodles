#![allow(dead_code)]

pub mod group;
mod parser;

pub use self::group::Group;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BaseModifications(Vec<Group>);

impl AsRef<[Group]> for BaseModifications {
    fn as_ref(&self) -> &[Group] {
        &self.0
    }
}

impl AsMut<Vec<Group>> for BaseModifications {
    fn as_mut(&mut self) -> &mut Vec<Group> {
        &mut self.0
    }
}

impl From<Vec<Group>> for BaseModifications {
    fn from(groups: Vec<Group>) -> Self {
        Self(groups)
    }
}

impl From<BaseModifications> for Vec<Group> {
    fn from(base_modifications: BaseModifications) -> Self {
        base_modifications.0
    }
}
