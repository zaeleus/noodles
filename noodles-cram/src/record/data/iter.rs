use std::{io, slice};

use noodles_sam::{self as sam, alignment::record::data::field::Tag};

use super::field::Value;

enum State<'r, 'c: 'r> {
    Fields(slice::Iter<'r, (Tag, Value<'c>)>),
    ReadGroup,
    Done,
}

pub(super) struct Iter<'r, 'c: 'r> {
    header: &'c sam::Header,
    state: State<'r, 'c>,
    read_group_id: Option<usize>,
}

impl<'r, 'c: 'r> Iter<'r, 'c> {
    pub(super) fn new(
        header: &'c sam::Header,
        fields: &'r [(Tag, Value<'c>)],
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

    fn size_hint(&self) -> (usize, Option<usize>) {
        match self.state {
            State::Fields(ref iter) => {
                let (lower, upper) = iter.size_hint();

                if self.read_group_id.is_none() {
                    (lower, upper)
                } else {
                    (lower + 1, upper.map(|n| n + 1))
                }
            }
            State::ReadGroup => {
                if self.read_group_id.is_none() {
                    (0, Some(0))
                } else {
                    (1, Some(1))
                }
            }
            State::Done => (0, Some(0)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_size_hint() {
        let header = sam::Header::default();
        assert_eq!(Iter::new(&header, &[], None).size_hint(), (0, Some(0)));
        assert_eq!(Iter::new(&header, &[], Some(0)).size_hint(), (1, Some(1)));
    }
}
