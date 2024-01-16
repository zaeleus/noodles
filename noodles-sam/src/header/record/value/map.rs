//! SAM header record map value.

pub mod builder;
pub mod header;
pub mod program;
pub mod read_group;
pub mod reference_sequence;
pub mod tag;

pub use self::{
    builder::Builder, header::Header, program::Program, read_group::ReadGroup,
    reference_sequence::ReferenceSequence, tag::Tag,
};

use indexmap::IndexMap;

pub(crate) type OtherFields<S> = IndexMap<tag::Other<S>, Vec<u8>>;

/// An inner SAM header record map value.
pub trait Inner: Sized {
    /// The standard tag type.
    type StandardTag: tag::Standard;

    /// The builder type.
    type Builder: builder::Inner<Self>;
}

/// A SAM header record map value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Map<I>
where
    I: Inner,
{
    pub(crate) inner: I,
    pub(crate) other_fields: OtherFields<I::StandardTag>,
}

impl<I> Map<I>
where
    I: Inner,
{
    /// Creates a SAM header record map value.
    pub fn builder() -> Builder<I> {
        Builder::default()
    }

    /// Returns the nonstandard fields in the map.
    pub fn other_fields(&self) -> &OtherFields<I::StandardTag> {
        &self.other_fields
    }

    /// Returns a mutable reference to the nonstandard fields in the map.
    ///
    /// # Example
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::{tag, Header}, Map};
    /// let mut map = Map::<Header>::new(Default::default());
    /// let nd = tag::Other::try_from([b'n', b'd'])?;
    /// map.other_fields_mut().insert(nd, Vec::from("noodles"));
    /// # Ok::<_, tag::ParseError>(())
    /// ```
    pub fn other_fields_mut(&mut self) -> &mut OtherFields<I::StandardTag> {
        &mut self.other_fields
    }
}

impl<I> Default for Map<I>
where
    I: Inner + Default,
{
    fn default() -> Self {
        Self {
            inner: I::default(),
            other_fields: OtherFields::new(),
        }
    }
}
