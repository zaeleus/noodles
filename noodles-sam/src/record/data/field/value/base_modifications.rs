//! SAM record data field value for base modifications.

pub mod group;
mod parser;

pub use self::group::Group;

use crate::alignment::record_buf::Sequence;

/// Base modifications.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BaseModifications(Vec<Group>);

impl BaseModifications {
    /// Parses base modifications from a string.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     alignment::record_buf::Sequence,
    ///     record::data::field::value::{
    ///         base_modifications::{
    ///             group::{modification, Strand, UnmodifiedBase},
    ///             Group,
    ///         },
    ///         BaseModifications,
    ///     },
    /// };
    ///
    /// let is_reverse_complemented = false;
    /// let sequence = Sequence::from(b"CACCCGATGACCGGCT".to_vec());
    /// let base_modifications = BaseModifications::parse(
    ///     "C+m,1,3,0;",
    ///     is_reverse_complemented,
    ///     &sequence,
    /// )?;
    ///
    /// assert_eq!(base_modifications, BaseModifications::from(vec![
    ///     Group::new(
    ///         UnmodifiedBase::C,
    ///         Strand::Forward,
    ///         vec![modification::FIVE_METHYLCYTOSINE],
    ///         None,
    ///         vec![2, 11, 14],
    ///     )
    /// ]));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn parse(
        s: &str,
        is_reverse_complemented: bool,
        sequence: &Sequence,
    ) -> Result<Self, parser::ParseError> {
        parser::parse(s, is_reverse_complemented, sequence)
    }
}

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
