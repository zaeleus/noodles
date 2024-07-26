use bstr::{BStr, BString};

/// Feature record other fields.
#[derive(Clone, Default, Debug, Eq, PartialEq)]
pub struct OtherFields(Vec<Option<BString>>);

impl AsRef<[Option<BString>]> for OtherFields {
    fn as_ref(&self) -> &[Option<BString>] {
        &self.0
    }
}

impl AsMut<Vec<Option<BString>>> for OtherFields {
    fn as_mut(&mut self) -> &mut Vec<Option<BString>> {
        &mut self.0
    }
}

impl crate::feature::record::OtherFields for OtherFields {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = &BStr> + '_> {
        todo!()
    }
}
