use std::io;

use bstr::{BStr, BString, ByteVec};
use indexmap::{IndexMap, IndexSet};

use super::record::value::map::{program::tag, Map, Program};

type Inner = IndexMap<BString, Map<Program>>;

#[allow(clippy::tabs_in_doc_comments)]
/// SAM header programs.
///
/// SAM header programs are header records that form program chains. A program chain is a directed
/// acyclic graph (DAG) starting at a root program and ending at a leaf program. Program edges are
/// directed forward from a parent program to its linked program using the previous program ID
/// (`PP`) field in the child program.
///
/// For example, take the following program records:
///
/// ```text
/// @PG	ID:pg0
/// @PG	ID:pg1	PP:pg0
/// ```
///
/// This creates the program chain `pg0 -> pg1`.
///
/// There can exist more than one program chain in the SAM header programs, e.g.,
///
/// ```text
/// @PG	ID:pg0
/// @PG	ID:pg1
/// @PG	ID:pg2
/// @PG	ID:pg3	PP:pg0
/// @PG	ID:pg4	PP:pg1
/// @PG	ID:pg5	PP:pg4
/// ```
///
/// This creates the program chains `pg0 -> pg3`, `pg1 -> pg4 -> pg5`, and `pg2`.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Programs(Inner);

impl Programs {
    /// Adds a program.
    ///
    /// If the program is the first program in the graph, this is similar to calling
    /// `IndexMap::insert` on the inner graph. If no previous program is set, this attaches the
    /// program to all program chains using leaf programs as the given program's previous program.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     self as sam,
    ///     header::record::value::{map::{program::tag, Program}, Map},
    /// };
    ///
    /// let mut header = sam::Header::default();
    /// let programs = header.programs_mut();
    ///
    /// programs.add("pg0", Map::default())?;
    /// programs.add("pg1", Map::default())?;
    ///
    /// let expected = sam::Header::builder()
    ///     .add_program("pg0", Map::default())
    ///     .add_program("pg1", Map::builder().insert(tag::PREVIOUS_PROGRAM_ID, "pg0").build()?)
    ///     .build();
    ///
    /// assert_eq!(programs, expected.programs());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn add<P>(&mut self, id_prefix: P, map: Map<Program>) -> io::Result<()>
    where
        P: Into<BString>,
    {
        const SEPARATOR: u8 = b'-';

        let id_prefix = id_prefix.into();

        if self.0.is_empty() {
            self.0.insert(id_prefix, map);
            return Ok(());
        }

        let previous_program_ids: Vec<BString> = self.leaves()?.map(|(id, _)| id.into()).collect();
        let contains_prefix_id = self.0.contains_key(&id_prefix);

        for (i, previous_program_id) in previous_program_ids.into_iter().enumerate() {
            let mut id = id_prefix.clone();

            if i > 0 || contains_prefix_id {
                id.push_byte(SEPARATOR);
                id.push_str(&previous_program_id);

                if self.0.contains_key(&id) {
                    return Err(io::Error::new(io::ErrorKind::InvalidInput, "duplicate ID"));
                }
            }

            let mut map = map.clone();

            map.other_fields_mut()
                .entry(tag::PREVIOUS_PROGRAM_ID)
                .or_insert(previous_program_id);

            self.0.insert(id, map);
        }

        Ok(())
    }

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
    fn test_add() -> Result<(), Box<dyn std::error::Error>> {
        let mut programs = Programs::default();
        assert!(programs.as_ref().is_empty());

        programs.add("pg0", Map::default())?;
        let expected = Programs(
            [(BString::from("pg0"), Map::default())]
                .into_iter()
                .collect(),
        );
        assert_eq!(programs, expected);

        programs.add("pg1", Map::default())?;
        let expected = Programs(
            [
                (BString::from("pg0"), Map::default()),
                (
                    BString::from("pg1"),
                    Map::builder()
                        .insert(tag::PREVIOUS_PROGRAM_ID, "pg0")
                        .build()?,
                ),
            ]
            .into_iter()
            .collect(),
        );
        assert_eq!(programs, expected);

        programs
            .as_mut()
            .insert(BString::from("pg2"), Map::default());
        programs.add("pg3", Map::default())?;
        let expected = Programs(
            [
                (BString::from("pg0"), Map::default()),
                (
                    BString::from("pg1"),
                    Map::builder()
                        .insert(tag::PREVIOUS_PROGRAM_ID, "pg0")
                        .build()?,
                ),
                (BString::from("pg2"), Map::default()),
                (
                    BString::from("pg3"),
                    Map::builder()
                        .insert(tag::PREVIOUS_PROGRAM_ID, "pg2")
                        .build()?,
                ),
                (
                    BString::from("pg3-pg1"),
                    Map::builder()
                        .insert(tag::PREVIOUS_PROGRAM_ID, "pg1")
                        .build()?,
                ),
            ]
            .into_iter()
            .collect(),
        );
        assert_eq!(programs, expected);

        Ok(())
    }

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
