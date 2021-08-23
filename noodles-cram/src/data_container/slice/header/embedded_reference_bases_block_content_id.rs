use std::ops::Deref;

const UNSET: i32 = -1;
const MIN: i32 = 0;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct EmbeddedReferenceBasesBlockContentId(Option<i32>);

impl Deref for EmbeddedReferenceBasesBlockContentId {
    type Target = Option<i32>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<i32> for EmbeddedReferenceBasesBlockContentId {
    fn from(n: i32) -> Self {
        if n < MIN {
            Self(None)
        } else {
            Self(Some(n))
        }
    }
}

impl From<EmbeddedReferenceBasesBlockContentId> for i32 {
    fn from(block_content_id: EmbeddedReferenceBasesBlockContentId) -> Self {
        match *block_content_id {
            Some(id) => id,
            None => UNSET,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_i32_for_embedded_reference_bases_block_content_id() {
        assert_eq!(*EmbeddedReferenceBasesBlockContentId::from(-1), None);
        assert_eq!(*EmbeddedReferenceBasesBlockContentId::from(13), Some(13));
    }

    #[test]
    fn test_from_embedded_reference_bases_block_content_id_for_i32() {
        assert_eq!(
            i32::from(EmbeddedReferenceBasesBlockContentId::from(-1)),
            -1
        );
        assert_eq!(
            i32::from(EmbeddedReferenceBasesBlockContentId::from(13)),
            13
        );
    }
}
