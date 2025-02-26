use std::{io, marker::PhantomData, mem};

use noodles_sam as sam;

/// BAM record data field array values.
#[derive(Debug, PartialEq)]
pub struct Values<'a, N> {
    src: &'a [u8],
    _marker: PhantomData<N>,
}

impl<'a, N> Values<'a, N> {
    pub(crate) fn new(src: &'a [u8]) -> Self {
        Self {
            src,
            _marker: PhantomData,
        }
    }
}

impl Values<'_, i8> {
    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<i8>> + '_ {
        self.src.iter().map(|&b| Ok(b as i8))
    }
}

impl<'a> sam::alignment::record::data::field::value::array::Values<'a, i8> for Values<'a, i8> {
    fn len(&self) -> usize {
        self.src.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<i8>> + '_> {
        Box::new(self.iter())
    }
}

impl Values<'_, u8> {
    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<u8>> + '_ {
        self.src.iter().copied().map(Ok)
    }
}

impl<'a> sam::alignment::record::data::field::value::array::Values<'a, u8> for Values<'a, u8> {
    fn len(&self) -> usize {
        self.src.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u8>> + '_> {
        Box::new(self.iter())
    }
}

impl Values<'_, i16> {
    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<i16>> + '_ {
        self.src.chunks(mem::size_of::<i16>()).map(|chunk| {
            let buf = chunk
                .try_into()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(i16::from_le_bytes(buf))
        })
    }
}

impl<'a> sam::alignment::record::data::field::value::array::Values<'a, i16> for Values<'a, i16> {
    fn len(&self) -> usize {
        self.src.len() / mem::size_of::<i16>()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<i16>> + '_> {
        Box::new(self.iter())
    }
}

impl Values<'_, u16> {
    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<u16>> + '_ {
        self.src.chunks(mem::size_of::<u16>()).map(|chunk| {
            let buf = chunk
                .try_into()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(u16::from_le_bytes(buf))
        })
    }
}

impl<'a> sam::alignment::record::data::field::value::array::Values<'a, u16> for Values<'a, u16> {
    fn len(&self) -> usize {
        self.src.len() / mem::size_of::<u16>()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u16>> + '_> {
        Box::new(self.iter())
    }
}

impl Values<'_, i32> {
    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<i32>> + '_ {
        self.src.chunks(mem::size_of::<i32>()).map(|chunk| {
            let buf = chunk
                .try_into()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(i32::from_le_bytes(buf))
        })
    }
}

impl<'a> sam::alignment::record::data::field::value::array::Values<'a, i32> for Values<'a, i32> {
    fn len(&self) -> usize {
        self.src.len() / mem::size_of::<i32>()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<i32>> + '_> {
        Box::new(self.iter())
    }
}

impl Values<'_, u32> {
    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<u32>> + '_ {
        self.src.chunks(mem::size_of::<u32>()).map(|chunk| {
            let buf = chunk
                .try_into()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(u32::from_le_bytes(buf))
        })
    }
}

impl<'a> sam::alignment::record::data::field::value::array::Values<'a, u32> for Values<'a, u32> {
    fn len(&self) -> usize {
        self.src.len() / mem::size_of::<u32>()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u32>> + '_> {
        Box::new(self.iter())
    }
}

impl Values<'_, f32> {
    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<f32>> + '_ {
        self.src.chunks(mem::size_of::<f32>()).map(|chunk| {
            let buf = chunk
                .try_into()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(f32::from_le_bytes(buf))
        })
    }
}

impl<'a> sam::alignment::record::data::field::value::array::Values<'a, f32> for Values<'a, f32> {
    fn len(&self) -> usize {
        self.src.len() / mem::size_of::<f32>()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<f32>> + '_> {
        Box::new(self.iter())
    }
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record::data::field::value::array::Values as _;

    use super::*;

    #[test]
    fn test_len() {
        assert_eq!(Values::<'_, i8>::new(&[]).len(), 0);
        assert_eq!(Values::<'_, i8>::new(&[0x00]).len(), 1);
        assert_eq!(Values::<'_, i8>::new(&[0x00, 0x00]).len(), 2);

        assert_eq!(Values::<'_, u8>::new(&[]).len(), 0);
        assert_eq!(Values::<'_, u8>::new(&[0x00]).len(), 1);
        assert_eq!(Values::<'_, u8>::new(&[0x00, 0x00]).len(), 2);

        assert_eq!(Values::<'_, i16>::new(&[]).len(), 0);
        assert_eq!(Values::<'_, i16>::new(&[0x00, 0x00]).len(), 1);
        assert_eq!(Values::<'_, i16>::new(&[0x00, 0x00, 0x00, 0x00]).len(), 2);

        assert_eq!(Values::<'_, u16>::new(&[]).len(), 0);
        assert_eq!(Values::<'_, u16>::new(&[0x00, 0x00]).len(), 1);
        assert_eq!(Values::<'_, u16>::new(&[0x00, 0x00, 0x00, 0x00]).len(), 2);

        assert_eq!(Values::<'_, i32>::new(&[]).len(), 0);
        assert_eq!(Values::<'_, i32>::new(&[0x00, 0x00, 0x00, 0x00]).len(), 1);
        assert_eq!(
            Values::<'_, i32>::new(&[0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]).len(),
            2
        );

        assert_eq!(Values::<'_, u32>::new(&[]).len(), 0);
        assert_eq!(Values::<'_, u32>::new(&[0x00, 0x00, 0x00, 0x00]).len(), 1);
        assert_eq!(
            Values::<'_, u32>::new(&[0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]).len(),
            2
        );

        assert_eq!(Values::<'_, f32>::new(&[]).len(), 0);
        assert_eq!(Values::<'_, f32>::new(&[0x00, 0x00, 0x00, 0x00]).len(), 1);
        assert_eq!(
            Values::<'_, f32>::new(&[0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]).len(),
            2
        );
    }

    #[test]
    fn test_iter() -> io::Result<()> {
        assert!(Values::<'_, i8>::new(&[]).iter().next().is_none());
        assert_eq!(
            Values::<'_, i8>::new(&[0x05])
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [5]
        );
        assert_eq!(
            Values::<'_, i8>::new(&[0xf8, 0x0d])
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [-8, 13]
        );

        assert!(Values::<'_, u8>::new(&[]).iter().next().is_none());
        assert_eq!(
            Values::<'_, u8>::new(&[0x05])
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [5]
        );
        assert_eq!(
            Values::<'_, u8>::new(&[0x08, 0x0d])
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [8, 13]
        );

        assert!(Values::<'_, i16>::new(&[]).iter().next().is_none());
        assert_eq!(
            Values::<'_, i16>::new(&[0x05, 0x00])
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [5]
        );
        assert_eq!(
            Values::<'_, i16>::new(&[0xf8, 0xff, 0x0d, 0x00])
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [-8, 13]
        );

        assert!(Values::<'_, u16>::new(&[]).iter().next().is_none());
        assert_eq!(
            Values::<'_, u16>::new(&[0x05, 0x00])
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [5]
        );
        assert_eq!(
            Values::<'_, u16>::new(&[0x08, 0x00, 0x0d, 0x00])
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [8, 13]
        );

        assert!(Values::<'_, i32>::new(&[]).iter().next().is_none());
        assert_eq!(
            Values::<'_, i32>::new(&[0x05, 0x00, 0x00, 0x00])
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [5]
        );
        assert_eq!(
            Values::<'_, i32>::new(&[0xf8, 0xff, 0xff, 0xff, 0x0d, 0x00, 0x00, 0x00])
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [-8, 13]
        );

        assert!(Values::<'_, u32>::new(&[]).iter().next().is_none());
        assert_eq!(
            Values::<'_, u32>::new(&[0x05, 0x00, 0x00, 0x00])
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [5]
        );
        assert_eq!(
            Values::<'_, u32>::new(&[0x08, 0x00, 0x00, 0x00, 0x0d, 0x00, 0x00, 0x00])
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [8, 13]
        );

        assert!(Values::<'_, f32>::new(&[]).iter().next().is_none());
        assert_eq!(
            Values::<'_, f32>::new(&[0x00, 0x00, 0x00, 0x00])
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [0.0]
        );
        assert_eq!(
            Values::<'_, f32>::new(&[0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00])
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [0.0, 0.0]
        );

        Ok(())
    }
}
