mod iter;

use std::io;

use noodles_sam::{
    self as sam,
    alignment::{record::data::field::Tag, record_buf::data::field::Value as ValueBuf},
};

use self::iter::Iter;

pub struct Data<'r, 'c: 'r> {
    header: &'c sam::Header,
    fields: &'r [(Tag, ValueBuf)],
    read_group_id: Option<usize>,
}

impl<'r, 'c: 'r> Data<'r, 'c> {
    pub(super) fn new(
        header: &'c sam::Header,
        fields: &'r [(Tag, ValueBuf)],
        read_group_id: Option<usize>,
    ) -> Self {
        Self {
            header,
            fields,
            read_group_id,
        }
    }
}

impl<'r, 'c: 'r> sam::alignment::record::Data for Data<'r, 'c> {
    fn is_empty(&self) -> bool {
        self.fields.is_empty()
    }

    fn get(&self, tag: &Tag) -> Option<io::Result<sam::alignment::record::data::field::Value<'_>>> {
        if *tag == Tag::READ_GROUP {
            return self.read_group_id.map(|id| {
                self.header
                    .read_groups()
                    .get_index(id)
                    .map(|(name, _)| {
                        sam::alignment::record::data::field::Value::String(name.as_ref())
                    })
                    .ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "invalid read group ID")
                    })
            });
        }

        for result in self.iter() {
            match result {
                Ok((t, value)) => {
                    if t == *tag {
                        return Some(Ok(value));
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        None
    }

    fn iter(
        &self,
    ) -> Box<
        dyn Iterator<Item = io::Result<(Tag, sam::alignment::record::data::field::Value<'_>)>> + '_,
    > {
        Box::new(Iter::new(self.header, self.fields, self.read_group_id))
    }
}
