use super::{
    record::value::{
        map::{self, Program, ReadGroup, ReferenceSequence},
        Map,
    },
    Header, Programs, ReadGroups, ReferenceSequences,
};

/// A SAM header builder.
#[derive(Debug, Default)]
pub struct Builder {
    header: Option<Map<map::Header>>,
    reference_sequences: ReferenceSequences,
    read_groups: ReadGroups,
    programs: Programs,
    comments: Vec<Vec<u8>>,
}

impl Builder {
    /// Sets a SAM header header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let header = sam::Header::builder()
    ///     .set_header(Default::default())
    ///     .build();
    ///
    /// assert!(header.header().is_some());
    /// ```
    pub fn set_header(mut self, header: Map<map::Header>) -> Self {
        self.header = Some(header);
        self
    }

    /// Sets the reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::{
    ///     self as sam,
    ///     header::record::value::{
    ///         map::{reference_sequence::Name, ReferenceSequence},
    ///         Map,
    ///     }
    /// };
    ///
    /// let reference_sequences = [(
    ///     Vec::from("sq0"),
    ///     Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?),
    /// )]
    /// .into_iter()
    /// .collect();
    ///
    /// let header = sam::Header::builder()
    ///     .set_reference_sequences(reference_sequences)
    ///     .build();
    ///
    /// let reference_sequences = header.reference_sequences();
    /// assert_eq!(reference_sequences.len(), 1);
    /// assert!(reference_sequences.contains_key(&b"sq0"[..]));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn set_reference_sequences(mut self, reference_sequences: ReferenceSequences) -> Self {
        self.reference_sequences = reference_sequences;
        self
    }

    /// Adds a reference sequence to the SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZeroUsize;
    ///
    /// use noodles_sam::{
    ///     self as sam,
    ///     header::record::value::{map::ReferenceSequence, Map},
    /// };
    ///
    /// let header = sam::Header::builder()
    ///     .add_reference_sequence(
    ///         "sq0",
    ///         Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?),
    ///     )
    ///     .build();
    ///
    /// let reference_sequences = header.reference_sequences();
    /// assert_eq!(reference_sequences.len(), 1);
    /// assert!(reference_sequences.contains_key(&b"sq0"[..]));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn add_reference_sequence<N>(
        mut self,
        name: N,
        reference_sequence: Map<ReferenceSequence>,
    ) -> Self
    where
        N: Into<Vec<u8>>,
    {
        self.reference_sequences
            .insert(name.into(), reference_sequence);

        self
    }

    /// Adds a read group to the SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{
    ///     self as sam,
    ///     header::record::value::{map::ReadGroup, Map},
    /// };
    ///
    /// let header = sam::Header::builder()
    ///     .add_read_group("rg0", Map::<ReadGroup>::default())
    ///     .build();
    ///
    /// let read_groups = header.read_groups();
    /// assert_eq!(read_groups.len(), 1);
    /// assert!(read_groups.contains_key(&b"rg0"[..]));
    /// ```
    pub fn add_read_group<I>(mut self, id: I, map: Map<ReadGroup>) -> Self
    where
        I: Into<Vec<u8>>,
    {
        self.read_groups.insert(id.into(), map);
        self
    }

    /// Adds a program to the SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::record::value::{map::Program, Map}};
    ///
    /// let header = sam::Header::builder()
    ///     .add_program("noodles-sam", Map::<Program>::default())
    ///     .build();
    ///
    /// let programs = header.programs();
    /// assert_eq!(programs.len(), 1);
    /// assert!(programs.contains_key(&b"noodles-sam"[..]));
    /// ```
    pub fn add_program<I>(mut self, id: I, map: Map<Program>) -> Self
    where
        I: Into<Vec<u8>>,
    {
        self.programs.insert(id.into(), map);
        self
    }

    /// Adds a comment to the SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let header = sam::Header::builder().add_comment("noodles-sam").build();
    /// let comments = header.comments();
    /// assert_eq!(comments.len(), 1);
    /// assert_eq!(&comments[0], b"noodles-sam");
    /// ```
    pub fn add_comment<C>(mut self, comment: C) -> Self
    where
        C: Into<Vec<u8>>,
    {
        self.comments.push(comment.into());
        self
    }

    /// Builds a SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let header = sam::Header::builder().build();
    /// assert!(header.is_empty());
    /// ```
    pub fn build(self) -> Header {
        Header {
            header: self.header,
            reference_sequences: self.reference_sequences,
            read_groups: self.read_groups,
            programs: self.programs,
            comments: self.comments,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let header = Builder::default();

        assert!(header.header.is_none());
        assert!(header.reference_sequences.is_empty());
        assert!(header.read_groups.is_empty());
        assert!(header.programs.is_empty());
        assert!(header.comments.is_empty());
    }

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZeroUsize;

        let header = Builder::default()
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(8)?),
            )
            .add_reference_sequence(
                "sq1",
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?),
            )
            .add_reference_sequence(
                "sq2",
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(21)?),
            )
            .add_read_group("rg0", Map::<ReadGroup>::default())
            .add_read_group("rg1", Map::<ReadGroup>::default())
            .add_program("pg0", Map::<Program>::default())
            .add_comment("written by noodles-sam")
            .build();

        let reference_sequences = header.reference_sequences();
        assert_eq!(reference_sequences.len(), 3);
        assert!(reference_sequences.contains_key(&b"sq0"[..]));
        assert!(reference_sequences.contains_key(&b"sq1"[..]));
        assert!(reference_sequences.contains_key(&b"sq2"[..]));

        assert_eq!(header.read_groups().len(), 2);

        assert_eq!(header.programs().len(), 1);

        let comments = header.comments();
        assert_eq!(comments.len(), 1);
        assert_eq!(&comments[0], b"written by noodles-sam");

        Ok(())
    }
}
