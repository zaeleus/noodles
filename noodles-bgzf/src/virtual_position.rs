#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct VirtualPosition(u64);

impl VirtualPosition {
    pub fn compressed(self) -> u64 {
        self.0 >> 16
    }

    pub fn uncompressed(self) -> u64 {
        self.0 & 0xffff
    }
}

impl From<u64> for VirtualPosition {
    fn from(pos: u64) -> Self {
        Self(pos)
    }
}

impl From<(u64, u64)> for VirtualPosition {
    fn from(pos: (u64, u64)) -> Self {
        Self(pos.0 << 16 | pos.1)
    }
}

impl From<VirtualPosition> for u64 {
    fn from(pos: VirtualPosition) -> Self {
        pos.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_u64_for_virtual_position() {
        let pos = VirtualPosition::from(88384945211);
        assert_eq!(pos.compressed(), 1348647);
        assert_eq!(pos.uncompressed(), 15419);

        let pos = VirtualPosition::from(188049630896);
        assert_eq!(pos.compressed(), 2869409);
        assert_eq!(pos.uncompressed(), 42672);

        let pos = VirtualPosition::from(26155658182977);
        assert_eq!(pos.compressed(), 399103671);
        assert_eq!(pos.uncompressed(), 321);
    }

    #[test]
    fn test_from_u64_u64_tuple_for_virtual_position() {
        assert_eq!(
            VirtualPosition::from((1348647, 15419)),
            VirtualPosition::from(88384945211)
        );

        assert_eq!(
            VirtualPosition::from((2869409, 42672)),
            VirtualPosition::from(188049630896)
        );

        assert_eq!(
            VirtualPosition::from((399103671, 321)),
            VirtualPosition::from(26155658182977)
        );
    }

    #[test]
    fn test_from_virtual_position_for_u64() {
        assert_eq!(u64::from(VirtualPosition::from(88384945211)), 88384945211);
        assert_eq!(u64::from(VirtualPosition::from(188049630896)), 188049630896);
        assert_eq!(
            u64::from(VirtualPosition::from(26155658182977)),
            26155658182977
        );
    }
}
