use std::io::{self, Write};

pub fn encode(src: &[u8], n: usize) -> io::Result<(Vec<Vec<u32>>, Vec<u8>)> {
    let contexts = build_contexts(src, n);
    let freq = normalize_contexts(contexts);
    let _cfreq = build_cumulative_contexts(&freq);

    Ok((freq, Vec::new()))
}

pub fn write_frequencies<W>(_writer: &mut W, _frequencies: &[Vec<u32>]) -> io::Result<()>
where
    W: Write,
{
    todo!()
}

fn build_contexts(src: &[u8], n: usize) -> Vec<Vec<u32>> {
    let mut frequencies = vec![vec![0; 256]; 256];

    let quarter = src.len() / n;

    for i in 0..n {
        let sym = usize::from(src[i * quarter]);
        frequencies[0][sym] += 1;
    }

    for window in src.windows(2) {
        let sym_0 = usize::from(window[0]);
        let sym_1 = usize::from(window[1]);
        frequencies[sym_0][sym_1] += 1;
    }

    frequencies
}

fn normalize_contexts(contexts: Vec<Vec<u32>>) -> Vec<Vec<u32>> {
    use super::normalize_frequencies_nx16_0;

    contexts
        .into_iter()
        .map(|frequencies| normalize_frequencies_nx16_0(&frequencies))
        .collect()
}

fn build_cumulative_contexts(contexts: &[Vec<u32>]) -> Vec<Vec<u32>> {
    use super::build_cumulative_frequencies;

    contexts
        .iter()
        .map(|frequencies| build_cumulative_frequencies(frequencies))
        .collect()
}
