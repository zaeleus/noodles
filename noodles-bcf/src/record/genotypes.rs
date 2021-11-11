#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Genotypes(Vec<u8>);

impl AsRef<[u8]> for Genotypes {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<Vec<u8>> for Genotypes {
    fn as_mut(&mut self) -> &mut Vec<u8> {
        &mut self.0
    }
}
