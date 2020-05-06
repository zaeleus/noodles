use super::{header, Header, Program, ReadGroup, ReferenceSequence, ReferenceSequences};

#[derive(Debug, Default)]
pub struct Builder {
    reference_sequences: ReferenceSequences,
    read_groups: Vec<ReadGroup>,
    programs: Vec<Program>,
    comments: Vec<String>,
}

impl Builder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_reference_sequence(mut self, reference_sequence: ReferenceSequence) -> Self {
        let name = reference_sequence.name().into();
        self.reference_sequences.insert(name, reference_sequence);
        self
    }

    pub fn add_read_group(mut self, read_group: ReadGroup) -> Self {
        self.read_groups.push(read_group);
        self
    }

    pub fn add_program(mut self, program: Program) -> Self {
        self.programs.push(program);
        self
    }

    pub fn add_comment<S>(mut self, comment: S) -> Self
    where
        S: Into<String>,
    {
        self.comments.push(comment.into());
        self
    }

    pub fn build(self) -> Header {
        Header {
            header: header::Header::default(),
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
        let header = Builder::default().build();

        assert_eq!(header.header().version(), "1.6");
        assert!(header.reference_sequences().is_empty());
        assert!(header.read_groups().is_empty());
        assert!(header.programs().is_empty());
        assert!(header.comments().is_empty());
    }

    #[test]
    fn test_build() {
        let header = Builder::new()
            .add_reference_sequence(ReferenceSequence::new(String::from("sq0"), 8))
            .add_reference_sequence(ReferenceSequence::new(String::from("sq1"), 13))
            .add_reference_sequence(ReferenceSequence::new(String::from("sq2"), 21))
            .add_read_group(ReadGroup::new(String::from("rg0")))
            .add_read_group(ReadGroup::new(String::from("rg1")))
            .add_program(Program::new(String::from("noodles-sam")))
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
