use super::{build_cumulative_frequencies, normalize_frequencies};

#[allow(dead_code)]
pub fn build_contexts(data: &[u8], bin_count: usize) -> Vec<Vec<u32>> {
    let mut frequencies = vec![vec![0; bin_count]; bin_count];

    let quarter = data.len() / 4;

    for i in 0..4 {
        frequencies[0][data[i * quarter] as usize] += 1;
    }

    for window in data.windows(2) {
        let sym_0 = window[0] as usize;
        let sym_1 = window[1] as usize;
        frequencies[sym_0][sym_1] += 1;
    }

    frequencies
}

#[allow(dead_code)]
pub fn normalize_contexts(contexts: &[Vec<u32>]) -> Vec<Vec<u32>> {
    contexts
        .iter()
        .map(|frequencies| normalize_frequencies(frequencies))
        .collect()
}

#[allow(dead_code)]
pub fn build_cumulative_contexts(contexts: &[Vec<u32>]) -> Vec<Vec<u32>> {
    contexts
        .iter()
        .map(|frequencies| build_cumulative_frequencies(frequencies))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_contexts() {
        let data = [1, 2, 3, 1, 2, 1, 2];
        let actual = build_contexts(&data, 4);
        let expected = [[0, 2, 1, 1], [0, 0, 3, 0], [0, 1, 0, 1], [0, 1, 0, 0]];
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_normalize_contexts() {
        let contexts = [
            vec![0, 2, 1, 1],
            vec![0, 0, 3, 0],
            vec![0, 1, 0, 1],
            vec![0, 1, 0, 0],
        ];
        let actual = normalize_contexts(&contexts);

        let expected = [
            [0, 2049, 1023, 1023],
            [0, 0, 4095, 0],
            [0, 2047, 0, 2048],
            [0, 4095, 0, 0],
        ];

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_build_cumulative_contexts() {
        let normalized_contexts = [
            vec![0, 2049, 1023, 1023],
            vec![0, 0, 4095, 0],
            vec![0, 2047, 0, 2048],
            vec![0, 4095, 0, 0],
        ];
        let actual = build_cumulative_contexts(&normalized_contexts);

        let expected = [
            [0, 0, 2049, 3072],
            [0, 0, 0, 4095],
            [0, 0, 2047, 2047],
            [0, 0, 4095, 4095],
        ];

        assert_eq!(actual, expected);
    }
}
