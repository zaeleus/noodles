use std::{io, iter::FusedIterator, slice};

use bstr::BStr;
use noodles_sam::{
    self as sam,
    alignment::{record::data::field::Tag, record_buf::data::field::Value as ValueBuf},
};

enum State<'r> {
    Fields(slice::Iter<'r, (Tag, ValueBuf)>),
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
        fields: &'r [(Tag, ValueBuf)],
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
                    Some((tag, value_buf)) => {
                        let sam_value: sam::alignment::record::data::field::Value<'_> =
                            value_buf.into();
                        return Some(Ok((*tag, sam_value)));
                    }
                    None => self.state = State::ReadGroup,
                },
                State::ReadGroup => {
                    self.state = State::Done;

                    if let Some(id) = self.read_group_id {
                        return Some(get_read_group_name(self.header, id).map(|name| {
                            (
                                Tag::READ_GROUP,
                                sam::alignment::record::data::field::Value::String(name),
                            )
                        }));
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

impl<'r, 'c: 'r> FusedIterator for Iter<'r, 'c> {}

fn get_read_group_name(header: &sam::Header, read_group_id: usize) -> io::Result<&BStr> {
    header
        .read_groups()
        .get_index(read_group_id)
        .map(|(name, _)| name.as_ref())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid read group ID"))
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
