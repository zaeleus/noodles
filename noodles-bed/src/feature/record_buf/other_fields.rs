use bstr::{BStr, BString};

/// Feature record other fields.
#[derive(Clone, Default, Debug, Eq, PartialEq)]
pub struct OtherFields(Vec<BString>);

impl AsRef<[BString]> for OtherFields {
    fn as_ref(&self) -> &[BString] {
        &self.0
    }
}

impl AsMut<Vec<BString>> for OtherFields {
    fn as_mut(&mut self) -> &mut Vec<BString> {
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
        Box::new(self.0.iter().map(|buf| buf.as_ref()))
    }
}

impl crate::feature::record::OtherFields for &OtherFields {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = &BStr> + '_> {
        Box::new(self.0.iter().map(|buf| buf.as_ref()))
    }
}
