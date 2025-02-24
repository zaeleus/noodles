/// The number of bytes of context used to compute frequencies.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Order {
    /// Order-0.
    Zero,
    /// Order-1.
    One,
}
