use super::{
    header, Header, Program, Programs, ReadGroup, ReadGroups, ReferenceSequence, ReferenceSequences,
};

/// A SAM header builder.
#[derive(Debug, Default)]
pub struct Builder {
    header: Option<header::Header>,
    reference_sequences: ReferenceSequences,
    read_groups: ReadGroups,
    programs: Programs,
    comments: Vec<String>,
}

impl Builder {
    /// Creates a new SAM header builder.
    ///
    /// Typically, [`Header::builder`] is used instead of calling this.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let builder = sam::Header::builder();
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets a SAM header header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    ///
    /// let header = sam::Header::builder()
    ///     .set_header(sam::header::header::Header::default())
    ///     .build();
    ///
    /// assert!(header.header().is_some());
    /// ```
    pub fn set_header(mut self, header: header::Header) -> Self {
        self.header = Some(header);
        self
    }

    /// Sets the reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::{ReferenceSequence, ReferenceSequences}};
    ///
    /// let reference_sequences: ReferenceSequences = vec![
    ///     (String::from("sq0"), ReferenceSequence::new("sq0", 13))
    /// ]
    /// .into_iter()
    /// .collect();
    ///
    /// let header = sam::Header::builder()
    ///     .set_reference_sequences(reference_sequences)
    ///     .build();
    ///
    /// let reference_sequences = header.reference_sequences();
    /// assert_eq!(reference_sequences.len(), 1);
    /// assert!(reference_sequences.contains_key("sq0"));
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
    /// use noodles_sam::{self as sam, header::ReferenceSequence};
    ///
    /// let header = sam::Header::builder()
    ///     .add_reference_sequence(ReferenceSequence::new("sq0", 13))
    ///     .build();
    ///
    /// let reference_sequences = header.reference_sequences();
    /// assert_eq!(reference_sequences.len(), 1);
    /// assert!(reference_sequences.contains_key("sq0"));
    /// ```
    pub fn add_reference_sequence(mut self, reference_sequence: ReferenceSequence) -> Self {
        let name = reference_sequence.name().into();
        self.reference_sequences.insert(name, reference_sequence);
        self
    }

    /// Adds a read group to the SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::ReadGroup};
    ///
    /// let header = sam::Header::builder()
    ///     .add_read_group(ReadGroup::new("rg0"))
    ///     .build();
    ///
    /// let read_groups = header.read_groups();
    /// assert_eq!(read_groups.len(), 1);
    /// assert!(read_groups.contains_key("rg0"));
    /// ```
    pub fn add_read_group(mut self, read_group: ReadGroup) -> Self {
        self.read_groups.insert(read_group.id().into(), read_group);
        self
    }

    /// Adds a program to the SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, header::Program};
    ///
    /// let header = sam::Header::builder()
    ///     .add_program(Program::new("noodles-sam"))
    ///     .build();
    ///
    /// let programs = header.programs();
    /// assert_eq!(programs.len(), 1);
    /// assert!(programs.contains_key("noodles-sam"));
    /// ```
    pub fn add_program(mut self, program: Program) -> Self {
        self.programs.insert(program.id().into(), program);
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
    /// assert_eq!(&comments[0], "noodles-sam");
    /// ```
    pub fn add_comment<S>(mut self, comment: S) -> Self
    where
        S: Into<String>,
    {
        self.comments.push(comment.into());
        self
    }

    /// Builds a SAM header.
    ///
    /// # Example
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
    fn test_build() {
        let header = Builder::new()
            .add_reference_sequence(ReferenceSequence::new("sq0", 8))
            .add_reference_sequence(ReferenceSequence::new("sq1", 13))
            .add_reference_sequence(ReferenceSequence::new("sq2", 21))
            .add_read_group(ReadGroup::new("rg0"))
            .add_read_group(ReadGroup::new("rg1"))
            .add_program(Program::new("noodles-sam"))
            .add_comment("written by noodles-sam")
            .build();

        let reference_sequences = header.reference_sequences();
        assert_eq!(reference_sequences.len(), 3);
        assert!(reference_sequences.contains_key("sq0"));
        assert!(reference_sequences.contains_key("sq1"));
        assert!(reference_sequences.contains_key("sq2"));

        assert_eq!(header.read_groups().len(), 2);

        assert_eq!(header.programs().len(), 1);

        let comments = header.comments();
        assert_eq!(comments.len(), 1);
        assert_eq!(&comments[0], "written by noodles-sam");
    }
}
