use bstr::{BStr, BString};
use indexmap::IndexMap;

use super::record::value::map::{program::tag, Map, Program};

type Inner = IndexMap<BString, Map<Program>>;

/// SAM header programs.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Programs(Inner);

impl Programs {
    /// Returns an iterator over root programs.
    ///
    /// A root program is a first program of a program chain.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     self as sam,
    ///     header::record::value::{map::{program::tag, Program}, Map},
    /// };
    ///
    /// let header = sam::Header::builder()
    ///     .add_program("pg0", Map::default())
    ///     .add_program("pg1", Map::builder().insert(tag::PREVIOUS_PROGRAM_ID, "pg0").build()?)
    ///     .add_program("pg2", Map::builder().insert(tag::PREVIOUS_PROGRAM_ID, "pg1").build()?)
    ///     .add_program("pg3", Map::default())
    ///     .build();
    ///
    /// let mut roots = header.programs().roots();
    /// assert_eq!(roots.next().map(|(id, _)| id.as_ref()), Some(&b"pg0"[..]));
    /// assert_eq!(roots.next().map(|(id, _)| id.as_ref()), Some(&b"pg3"[..]));
    /// assert!(roots.next().is_none());
    /// # Ok::<_, sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn roots(&self) -> impl Iterator<Item = (&BStr, &Map<Program>)> {
        self.0
            .iter()
            .filter(|(_, map)| !map.other_fields().contains_key(&tag::PREVIOUS_PROGRAM_ID))
            .map(|(id, map)| (id.as_ref(), map))
    }
}

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
