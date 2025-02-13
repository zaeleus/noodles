use std::{io, slice};

use noodles_sam::{
    self as sam,
    alignment::{record::data::field::Tag, record_buf::data::field::Value},
};

enum State<'r> {
    Fields(slice::Iter<'r, (Tag, Value)>),
    ReadGroup,
    Done,
}

pub(super) struct Iter<'r, 'c: 'r> {
    header: &'c sam::Header,
    state: State<'r>,
    read_group_id: Option<usize>,
}

impl<'r, 'c: 'r> Iter<'r, 'c> {
    pub(super) fn new(
        header: &'c sam::Header,
        fields: &'r [(Tag, Value)],
        read_group_id: Option<usize>,
    ) -> Self {
        Self {
            header,
            state: State::Fields(fields.iter()),
            read_group_id,
        }
    }
}

impl<'r, 'c: 'r> Iterator for Iter<'r, 'c> {
    type Item = io::Result<(Tag, sam::alignment::record::data::field::Value<'r>)>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.state {
                State::Fields(ref mut iter) => match iter.next() {
                    Some((tag, value)) => return Some(Ok((*tag, value.into()))),
                    None => self.state = State::ReadGroup,
                },
                State::ReadGroup => {
                    self.state = State::Done;

                    if let Some(id) = self.read_group_id {
                        return Some(
                            self.header
                                .read_groups()
                                .get_index(id)
                                .map(|(name, _)| {
                                    (
                                        Tag::READ_GROUP,
                                        sam::alignment::record::data::field::Value::String(
                                            name.as_ref(),
                                        ),
                                    )
                                })
                                .ok_or_else(|| {
                                    io::Error::new(
                                        io::ErrorKind::InvalidData,
                                        "invalid read group ID",
                                    )
                                }),
                        );
                    }
                }
                State::Done => return None,
            }
        }
    }
}
