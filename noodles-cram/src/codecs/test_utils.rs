//! Shared codec test helpers.
//!
//! `fixtures/` holds streams encoded by htscodecs 1.6.7 (`rans_compress` and `rans_compress_4x16`)
//! with their decoded bytes: `qual.bin` for `*-qual-*` and [`sparse`] for `*-mid-*`,
//! `*-big-*` and `*-wide-*`. `r4x16-*` files are named by their flags byte.

/// A 64-bit LCG, matching the fixture generator's `next()`.
pub(crate) struct Rng(pub(crate) u64);

impl Rng {
    pub(crate) fn next(&mut self, n: usize) -> usize {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as usize % n
    }
}

/// Mostly `A`, with about every `k`-th byte drawn from the `n` letters from `A`.
///
/// `sparse(8000, 1024, 8)` is the `mid` fixture input, `sparse(70001, 2048, 4)` the `big` one, and
/// `sparse(4099, 16, 32)` the `wide` one.
pub(crate) fn sparse(len: usize, k: usize, n: usize) -> Vec<u8> {
    let mut rng = Rng(1);

    (0..len)
        .map(|_| {
            if rng.next(k) == 0 {
                b'A' + rng.next(n) as u8
            } else {
                b'A'
            }
        })
        .collect()
}
