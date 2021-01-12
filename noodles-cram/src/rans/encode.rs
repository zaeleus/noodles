// Base `b`.
#[allow(dead_code)]
const BASE: usize = 256;

#[allow(dead_code)]
fn build_frequencies(data: &[u8], bin_count: usize) -> Vec<u32> {
    let mut frequencies = vec![0; bin_count];

    for &b in data {
        let i = usize::from(b);
        frequencies[i] += 1;
    }

    frequencies
}

#[allow(dead_code)]
fn normalize_frequencies(frequencies: &[u32]) -> Vec<u32> {
    const SCALE: u32 = 4095;

    let mut sum = 0;
    let mut max = 0;
    let mut max_index = 0;

    for (i, &f) in frequencies.iter().enumerate() {
        if f >= max {
            max = f;
            max_index = i;
        }

        sum += f;
    }

    let mut normalized_sum = 0;
    let mut normalized_frequencies = vec![0; frequencies.len()];

    for (i, &f) in frequencies.iter().enumerate() {
        let normalized_frequency = f * SCALE / sum;
        normalized_frequencies[i] = normalized_frequency;
        normalized_sum += normalized_frequency;
    }

    // Because the calculation of `normalized_frequency` uses integer division (truncation), it's
    // possible that the sum of all the normalized frequencies is smaller than the scale value. In
    // this case, the difference is added to the last max value.
    if normalized_sum < SCALE {
        normalized_frequencies[max_index] += SCALE - normalized_sum;
    }

    normalized_frequencies
}

#[allow(dead_code)]
fn build_cumulative_frequencies(frequencies: &[u32]) -> Vec<u32> {
    let mut cumulative_frequencies = vec![0; frequencies.len()];

    for i in 0..frequencies.len() - 1 {
        cumulative_frequencies[i + 1] = cumulative_frequencies[i] + frequencies[i];
    }

    cumulative_frequencies
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_frequencies() {
        let data = [0, 1, 1, 2, 2, 2];
        assert_eq!(build_frequencies(&data, 4), [1, 2, 3, 0]);
    }

    #[test]
    fn test_normalize_frequencies() {
        let frequencies = [1, 2, 3, 0];
        assert_eq!(normalize_frequencies(&frequencies), [682, 1365, 2048, 0]);
    }

    #[test]
    fn test_build_cumulative_frequencies() {
        let frequencies = [682, 1365, 2048, 0];

        assert_eq!(
            build_cumulative_frequencies(&frequencies),
            [0, 682, 2047, 4095]
        );
    }
}
