use super::Base;

/// A 2D histogram of reference-read base substitutions.
///
/// Each bin is the frequency of a substitution from the reference base to a read base. A same base
/// substitution is not a substitution; thus, the main diagonal will always be 0s.
///
/// ```text
///                         read base
///                     A   C   G   T   N
///                   +---+---+---+---+---+
///                 A | 0 |   |   |   |   |
///                   +---+---+---+---+---+
///                 C |   | 0 |   |   |   |
///                   +---+---+---+---+---+
/// reference base  G |   |   | 0 |   |   |
///                   +---+---+---+---+---+
///                 T |   |   |   | 0 |   |
///                   +---+---+---+---+---+
///                 N |   |   |   |   | 0 |
///                   +---+---+---+---+---+
/// ```
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Histogram {
    bins: [[u64; 5]; 5],
}

impl Histogram {
    pub fn new(bins: [[u64; 5]; 5]) -> Self {
        Self { bins }
    }

    pub fn get_bins(&self, reference_base: Base) -> &[u64] {
        let i = reference_base as usize;
        &self.bins[i]
    }

    pub fn get(&self, reference_base: Base, read_base: Base) -> u64 {
        let i = reference_base as usize;
        let j = read_base as usize;
        self.bins[i][j]
    }

    pub fn hit(&mut self, reference_base: Base, read_base: Base) {
        let i = reference_base as usize;
        let j = read_base as usize;
        self.bins[i][j] += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_bins() {
        let mut histogram = Histogram::default();
        histogram.hit(Base::A, Base::C);
        assert_eq!(histogram.get_bins(Base::A), [0, 1, 0, 0, 0]);
    }

    #[test]
    fn test_get() {
        let mut histogram = Histogram::default();
        assert_eq!(histogram.get(Base::A, Base::C), 0);
        histogram.hit(Base::A, Base::C);
        assert_eq!(histogram.get(Base::A, Base::C), 1);
    }

    #[test]
    fn test_hit() {
        let mut histogram = Histogram::default();

        histogram.hit(Base::A, Base::C);
        histogram.hit(Base::G, Base::T);
        histogram.hit(Base::T, Base::A);

        histogram.hit(Base::C, Base::G);
        histogram.hit(Base::C, Base::G);

        assert_eq!(
            histogram.bins,
            [
                [0, 1, 0, 0, 0],
                [0, 0, 2, 0, 0],
                [0, 0, 0, 1, 0],
                [1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
            ]
        );
    }
}
