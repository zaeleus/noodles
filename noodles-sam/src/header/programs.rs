use std::io;

use bstr::{BStr, BString};
use indexmap::{IndexMap, IndexSet};

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

    /// Returns an iterator over leaf programs.
    ///
    /// A leaf program is the last program of a program chain.
    ///
    /// # Errors
    ///
    /// This returns an `io::Error` if any program chain has a cycle.
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
    /// let mut leaves = header.programs().leaves()?;
    /// assert_eq!(leaves.next().map(|(id, _)| id.as_ref()), Some(&b"pg3"[..]));
    /// assert_eq!(leaves.next().map(|(id, _)| id.as_ref()), Some(&b"pg2"[..]));
    /// assert!(leaves.next().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn leaves(&self) -> io::Result<impl Iterator<Item = (&BStr, &Map<Program>)>> {
        let mut ids: IndexSet<_> = self.0.keys().collect();

        for (id, map) in &self.0 {
            if let Some(previous_program_id) = map.other_fields().get(&tag::PREVIOUS_PROGRAM_ID) {
                if has_cycle(&self.0, previous_program_id.as_ref(), id.as_ref()) {
                    return Err(io::Error::new(io::ErrorKind::InvalidData, "cycle detected"));
                }

                ids.swap_remove(previous_program_id);
            }
        }

        Ok(ids.into_iter().map(|id| {
            // SAFETY: `id` is guaranteed to be in the set of keys.
            self.0
                .get_key_value(id)
                .map(|(i, map)| (i.as_ref(), map))
                .unwrap()
        }))
    }
}

fn has_cycle<'a>(graph: &'a Inner, mut parent_id: &'a BStr, node_id: &'a BStr) -> bool {
    loop {
        // SAFETY: `parent_id` is guaranteed to be in the graph.
        let node = &graph[parent_id];

        if let Some(previous_program_id) = node.other_fields().get(&tag::PREVIOUS_PROGRAM_ID) {
            parent_id = previous_program_id.as_ref();

            if parent_id == node_id {
                return true;
            }
        } else {
            return false;
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Header;

    #[test]
    fn test_leaves() -> Result<(), Box<dyn std::error::Error>> {
        let header = Header::builder()
            .add_program("pg0", Map::default())
            .add_program("pg1", Map::default())
            .add_program("pg2", Map::default())
            .add_program(
                "pg3",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pg0")
                    .build()?,
            )
            .add_program(
                "pg4",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pg1")
                    .build()?,
            )
            .add_program(
                "pg5",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pg4")
                    .build()?,
            )
            .build();

        let mut leaves = header.programs().leaves()?;
        assert_eq!(leaves.next().map(|(id, _)| id.as_ref()), Some(&b"pg5"[..]));
        assert_eq!(leaves.next().map(|(id, _)| id.as_ref()), Some(&b"pg3"[..]));
        assert_eq!(leaves.next().map(|(id, _)| id.as_ref()), Some(&b"pg2"[..]));
        assert!(leaves.next().is_none());

        Ok(())
    }

    #[test]
    fn test_leaves_with_cycle() -> Result<(), crate::header::record::value::map::builder::BuildError>
    {
        let header = Header::builder()
            .add_program(
                "pg0",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pg1")
                    .build()?,
            )
            .add_program(
                "pg1",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pg0")
                    .build()?,
            )
            .build();

        assert!(matches!(
            header.programs().leaves(),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
