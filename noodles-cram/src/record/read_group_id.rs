// § 8.4 Compression header block (2020-06-22): "Special value '-1' stands for no group."
pub(crate) const MISSING: i32 = -1;

/// A CRAM record read group ID.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct ReadGroupId(usize);

impl From<usize> for ReadGroupId {
    fn from(n: usize) -> Self {
        Self(n)
    }
}

impl From<ReadGroupId> for usize {
    fn from(read_group_id: ReadGroupId) -> Self {
        read_group_id.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_usize_for_read_group_id() {
        assert_eq!(ReadGroupId::from(0), ReadGroupId(0));
        assert_eq!(ReadGroupId::from(13), ReadGroupId(13));
    }

    #[test]
    fn test_from_read_group_id_for_usize() {
        assert_eq!(usize::from(ReadGroupId(0)), 0);
        assert_eq!(usize::from(ReadGroupId(13)), 13);
    }
}
