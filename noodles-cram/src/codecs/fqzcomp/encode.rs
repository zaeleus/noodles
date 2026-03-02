use std::{
    io::{self, Write},
    num::NonZero,
};

use super::{
    Models,
    parameters::{self, parameter},
};
use crate::{
    codecs::aac::RangeCoder,
    io::writer::num::{write_u8, write_u16_le, write_uint7},
};

/// Encodes quality scores using the fqzcomp codec.
///
/// `records` is a slice of `(read_length, is_reverse_strand)` tuples.
/// `src` is the concatenated quality score data for all records.
pub fn encode(records: &[(usize, bool)], src: &[u8]) -> io::Result<Vec<u8>> {
    let mut dst = Vec::new();

    let len =
        u32::try_from(src.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_uint7(&mut dst, len)?;

    let has_reverse = records.iter().any(|(_, rev)| *rev);
    let lens: Vec<usize> = records.iter().map(|(l, _)| *l).collect();

    // Reverse quality scores for reverse-strand reads before encoding
    let mut working_src;
    let src = if has_reverse {
        working_src = src.to_vec();
        let mut offset = 0;
        for &(rec_len, is_rev) in records {
            if is_rev && rec_len > 1 {
                working_src[offset..offset + rec_len].reverse();
            }
            offset += rec_len;
        }
        &working_src
    } else {
        src
    };

    let parameters = build_parameters(&lens, src, has_reverse);

    // Build q_hist: maps original quality -> model symbol (identity or inverse_qmap).
    let q_hist: Vec<[u8; 256]> = parameters
        .params
        .iter()
        .map(|param| {
            if let Some(ref inv) = param.inverse_qmap {
                *inv
            } else {
                let mut m = [0u8; 256];
                for (i, v) in m.iter_mut().enumerate() {
                    *v = i as u8;
                }
                m
            }
        })
        .collect();

    fqz_encode_params(&mut dst, &parameters)?;

    let mut range_coder = RangeCoder::default();
    let mut models = Models::new(parameters.symbol_count, parameters.selector_count());

    let mut p: usize = 0;
    let mut rec_num: usize = 0;

    let mut x: usize = 0;
    let mut last: u32 = 0;
    let mut qlast: u32 = 0;
    let mut selector: u8 = 0;
    let mut delta: u32 = 0;
    let mut prev_q: u8 = 0;

    let mut i = 0;

    while i < src.len() {
        if p == 0 {
            // Encode selector (HAVE_S_TAB or MULTI_PARAM)
            if let Some(ref mut sel_model) = models.sel {
                selector = parameters
                    .record_selectors
                    .as_ref()
                    .map(|sels| sels[rec_num])
                    .unwrap_or(0);
                sel_model.encode(&mut dst, &mut range_coder, selector)?;

                if parameters.gflags.contains(parameters::Flags::HAVE_S_TAB) {
                    x = usize::from(parameters.s_tab[usize::from(selector)]);
                } else {
                    x = usize::from(selector);
                }
            } else {
                x = usize::from(parameters.s_tab[0]);
            }

            let param = &parameters.params[x];
            let rec_len = lens[rec_num];

            if !param.flags.is_fixed_length() || rec_num == 0 {
                encode_length(&mut dst, &mut range_coder, &mut models, rec_len)?;
            }

            // DO_REV: encode strand flag per record
            if parameters.gflags.contains(parameters::Flags::DO_REV) {
                let rev_byte = if records[rec_num].1 { 1u8 } else { 0u8 };
                models.rev.encode(&mut dst, &mut range_coder, rev_byte)?;
            }

            // DO_DEDUP: check if current record duplicates the previous data
            if param.flags.has_duplicates() {
                let is_dup =
                    rec_num > 0 && i >= rec_len && src[i..i + rec_len] == src[i - rec_len..i];

                let dup_byte = if is_dup { 1u8 } else { 0u8 };
                models.dup.encode(&mut dst, &mut range_coder, dup_byte)?;

                if is_dup {
                    i += rec_len;
                    rec_num += 1;
                    continue;
                }
            }

            p = rec_len;
            last = u32::from(param.context);
            qlast = 0;
            delta = 0;
            prev_q = 0;

            rec_num += 1;
        }

        let q = src[i];
        let qq = q_hist[x][usize::from(q)];
        models.qual[usize::from(last as u16)].encode(&mut dst, &mut range_coder, qq)?;

        let param = &parameters.params[x];

        qlast = (qlast << param.q_shift)
            .overflowing_add(u32::from(param.q_tab[usize::from(qq)]))
            .0;
        last = u32::from(param.context);
        last += (qlast & ((1 << param.q_bits) - 1)) << param.q_loc;

        if param.flags.has_positions_table() {
            last += u32::from(param.p_tab[p.min(1023)]) << param.p_loc;
        }

        if param.flags.has_deltas_table() {
            let d = delta.min(255) as usize;
            last += u32::from(param.d_tab[d]) << param.d_loc;

            if prev_q != qq {
                delta += 1;
            }

            prev_q = qq;
        }

        if param.flags.has_selector() {
            last += u32::from(selector) << param.s_loc;
        }

        last &= 0xffff;
        p -= 1;
        i += 1;
    }

    range_coder.range_encode_end(&mut dst)?;

    Ok(dst)
}

struct Parameters {
    pub gflags: parameters::Flags,
    pub max_sel: u8,
    pub s_tab: Vec<u8>,
    pub params: Vec<Parameter>,
    pub symbol_count: NonZero<usize>,
    pub record_selectors: Option<Vec<u8>>,
}

impl Parameters {
    fn selector_count(&self) -> Option<NonZero<usize>> {
        if self.gflags.contains(parameters::Flags::HAVE_S_TAB) {
            NonZero::new(usize::from(self.max_sel) + 1)
        } else if self.gflags.contains(parameters::Flags::MULTI_PARAM) {
            NonZero::new(self.params.len())
        } else {
            None
        }
    }
}

struct Parameter {
    pub context: u16,
    pub flags: parameter::Flags,

    pub symbol_count: NonZero<usize>,

    pub q_bits: u8,
    pub q_shift: u8,
    pub q_loc: u8,

    pub s_loc: u8,

    pub p_loc: u8,

    pub d_loc: u8,

    pub q_tab: Vec<u8>,
    pub p_tab: Vec<u8>,
    pub d_tab: Vec<u8>,
    pub quality_map: Option<Vec<u8>>,
    pub inverse_qmap: Option<[u8; 256]>,
}

fn build_parameters(lens: &[usize], src: &[u8], has_reverse: bool) -> Parameters {
    // Set global flags
    let mut gflags = parameters::Flags::empty();
    if has_reverse {
        gflags |= parameters::Flags::DO_REV;
    }

    // Try multi-param grouping
    if let Some((group_assignments, n_groups)) = assign_record_groups(lens, 10) {
        return build_multi_parameters(lens, src, gflags, &group_assignments, n_groups);
    }

    let (param, effective_max_symbol) = build_single_parameter(lens, src, false);

    let last_i = 0u8;
    let s_tab = vec![last_i; 256];

    let global_symbol_count = NonZero::new(usize::from(effective_max_symbol) + 1).unwrap();

    Parameters {
        gflags,
        max_sel: 0,
        s_tab,
        params: vec![param],
        symbol_count: global_symbol_count,
        record_selectors: None,
    }
}

fn build_single_parameter(lens: &[usize], src: &[u8], has_selector: bool) -> (Parameter, u8) {
    let mut max_symbol = u8::MIN;
    let mut seen = [false; 256];

    for &b in src {
        max_symbol = max_symbol.max(b);
        seen[usize::from(b)] = true;
    }
    let distinct: Vec<u8> = (0u8..=255).filter(|&v| seen[usize::from(v)]).collect();
    let distinct_count = distinct.len();

    let (quality_map, inverse_qmap, effective_max_symbol) =
        if distinct_count > 0 && distinct_count <= 16 {
            // qmap: model_symbol -> original_quality (sorted distinct values)
            let qmap = distinct;
            // inverse_qmap: original_quality -> model_symbol
            let mut inv = [0u8; 256];
            for (model_sym, &orig_q) in qmap.iter().enumerate() {
                inv[usize::from(orig_q)] = model_sym as u8;
            }
            let eff_max = (distinct_count - 1) as u8;
            (Some(qmap), Some(inv), eff_max)
        } else {
            (None, None, max_symbol)
        };

    // Detect duplicate records
    let mut dup_count: usize = 0;
    if lens.len() > 1 {
        let mut offset = 0;
        for window in lens.windows(2) {
            let prev_len = window[0];
            let cur_len = window[1];
            let prev_start = offset;
            offset += prev_len;
            let cur_start = offset;

            if prev_len == cur_len
                && cur_start + cur_len <= src.len()
                && src[prev_start..prev_start + prev_len] == src[cur_start..cur_start + cur_len]
            {
                dup_count += 1;
            }
        }
    }

    let dup_fraction = if lens.len() > 1 {
        dup_count as f64 / (lens.len() - 1) as f64
    } else {
        0.0
    };

    let q_shift = 5;
    let q_bits = if q_shift > 4 { 9 } else { 8 };
    let p_bits = 7;
    let p_shift = i32::from(lens.first().copied().unwrap_or(0) > 128);

    let mut p_tab = vec![0; 1024];

    for (i, p) in p_tab.iter_mut().enumerate() {
        *p = ((1 << p_bits) - 1).min(i >> p_shift) as u8;
    }

    let mut flags = parameter::Flags::HAVE_PTAB;

    if lens.len() > 1 && lens.windows(2).all(|w| w[0] == w[1]) {
        flags |= parameter::Flags::DO_LEN;
    }

    // Enable dedup when > 5% of records are duplicates
    if dup_fraction > 0.05 {
        flags |= parameter::Flags::DO_DEDUP;
    }

    if quality_map.is_some() {
        flags |= parameter::Flags::HAVE_QMAP;
    }

    if has_selector {
        flags |= parameter::Flags::DO_SEL;
    }

    // Build delta table when there's sufficient data for it to help
    let (d_tab, d_loc) = if src.len() > 256 {
        let d_bits = 3;
        let mut dtab = vec![0u8; 256];
        for (i, d) in dtab.iter_mut().enumerate() {
            *d = ((1 << d_bits) - 1).min(i) as u8;
        }
        flags |= parameter::Flags::HAVE_DTAB;
        (dtab, 15u8)
    } else {
        (Vec::new(), 15u8)
    };

    // Build quality context quantization table
    let q_tab = build_quality_table(effective_max_symbol, q_bits);
    let is_identity = q_tab.iter().enumerate().all(|(i, &v)| v == i as u8);
    if !is_identity {
        flags |= parameter::Flags::HAVE_QTAB;
    }

    let param = Parameter {
        context: 0,
        flags,
        symbol_count: NonZero::new(usize::from(effective_max_symbol) + 1).unwrap(),
        q_bits,
        q_shift,
        q_loc: 7,
        s_loc: 15,
        p_loc: 0,
        d_loc,
        q_tab,
        p_tab,
        d_tab,
        quality_map,
        inverse_qmap,
    };

    (param, effective_max_symbol)
}

/// Builds a quality context quantization table mapping quality values to context bins.
/// The table is non-decreasing (required by write_array/read_array RLE encoding).
fn build_quality_table(max_q: u8, q_bits: u8) -> Vec<u8> {
    let max_bin = ((1u16 << q_bits) - 1).min(255) as u8;

    if max_q == 0 {
        return (0..=u8::MAX).map(|_| 0u8).collect();
    }

    (0..=u8::MAX)
        .map(|i| {
            let bin =
                (u16::from(i) * u16::from(max_bin) / u16::from(max_q)).min(u16::from(max_bin));
            bin as u8
        })
        .collect()
}

/// Assigns each record to a parameter group by read length.
/// Returns None if multi-param isn't beneficial.
fn assign_record_groups(lens: &[usize], min_group_size: usize) -> Option<(Vec<u8>, usize)> {
    if lens.len() < min_group_size * 2 {
        return None;
    }

    // Find median read length
    let mut sorted_lens: Vec<usize> = lens.to_vec();
    sorted_lens.sort_unstable();
    let median = sorted_lens[sorted_lens.len() / 2];

    // Split into 2 groups: short (<= median) and long (> median)
    let assignments: Vec<u8> = lens
        .iter()
        .map(|&l| if l > median { 1 } else { 0 })
        .collect();

    let count0 = assignments.iter().filter(|&&a| a == 0).count();
    let count1 = assignments.iter().filter(|&&a| a == 1).count();

    if count0 < min_group_size || count1 < min_group_size {
        return None;
    }

    Some((assignments, 2))
}

fn build_multi_parameters(
    lens: &[usize],
    src: &[u8],
    mut gflags: parameters::Flags,
    group_assignments: &[u8],
    n_groups: usize,
) -> Parameters {
    gflags |= parameters::Flags::MULTI_PARAM;
    gflags |= parameters::Flags::HAVE_S_TAB;

    let max_sel = (n_groups - 1) as u8;

    // Build s_tab: [0, 1, ..., n_groups-1, n_groups-1, ..., n_groups-1]
    let mut s_tab = vec![max_sel; 256];
    for (i, entry) in s_tab.iter_mut().enumerate().take(n_groups) {
        *entry = i as u8;
    }

    // Split data by group
    let mut group_lens: Vec<Vec<usize>> = vec![Vec::new(); n_groups];
    let mut group_src: Vec<Vec<u8>> = vec![Vec::new(); n_groups];

    let mut offset = 0;
    for (rec_idx, &rec_len) in lens.iter().enumerate() {
        let g = usize::from(group_assignments[rec_idx]);
        group_lens[g].push(rec_len);
        group_src[g].extend_from_slice(&src[offset..offset + rec_len]);
        offset += rec_len;
    }

    // Build a parameter for each group
    let mut params = Vec::with_capacity(n_groups);
    let mut global_max_symbol: u8 = 0;

    for g in 0..n_groups {
        let (param, eff_max) = build_single_parameter(&group_lens[g], &group_src[g], true);
        global_max_symbol = global_max_symbol.max(eff_max);
        params.push(param);
    }

    let global_symbol_count = NonZero::new(usize::from(global_max_symbol) + 1).unwrap();

    Parameters {
        gflags,
        max_sel,
        s_tab,
        params,
        symbol_count: global_symbol_count,
        record_selectors: Some(group_assignments.to_vec()),
    }
}

fn fqz_encode_params<W>(writer: &mut W, parameters: &Parameters) -> io::Result<()>
where
    W: Write,
{
    const VERSION: u8 = 5;

    write_u8(writer, VERSION)?;

    let gflags = u8::from(parameters.gflags);
    write_u8(writer, gflags)?;

    if parameters.gflags.contains(parameters::Flags::MULTI_PARAM) {
        let n_param = u8::try_from(parameters.params.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        write_u8(writer, n_param)?;
    }

    if parameters.gflags.contains(parameters::Flags::HAVE_S_TAB) {
        write_u8(writer, parameters.max_sel)?;
        write_array(writer, &parameters.s_tab)?;
    }

    for param in &parameters.params {
        fqz_encode_single_param(writer, param)?;
    }

    Ok(())
}

fn fqz_encode_single_param<W>(writer: &mut W, parameter: &Parameter) -> io::Result<()>
where
    W: Write,
{
    let context = parameter.context;
    write_u16_le(writer, context)?;

    let pflags = u8::from(parameter.flags);
    write_u8(writer, pflags)?;

    let max_sym = (usize::from(parameter.symbol_count) - 1) as u8;
    write_u8(writer, max_sym)?;

    write_u8(writer, (parameter.q_bits << 4) | parameter.q_shift)?;
    write_u8(writer, (parameter.q_loc << 4) | parameter.s_loc)?;
    write_u8(writer, (parameter.p_loc << 4) | parameter.d_loc)?;

    if parameter.flags.has_quality_map()
        && let Some(ref qmap) = parameter.quality_map
    {
        writer.write_all(qmap)?;
    }

    if parameter.flags.has_qualities_table() {
        write_array(writer, &parameter.q_tab)?;
    }

    if parameter.flags.has_positions_table() {
        write_array(writer, &parameter.p_tab)?;
    }

    if parameter.flags.has_deltas_table() {
        write_array(writer, &parameter.d_tab)?;
    }

    Ok(())
}

fn write_array<W>(writer: &mut W, data: &[u8]) -> io::Result<()>
where
    W: Write,
{
    let mut rle1 = Vec::new();

    let mut i = 0;
    let mut j = 0;

    while j < data.len() {
        let start = j;

        while j < data.len() && usize::from(data[j]) == i {
            j += 1;
        }

        let mut len = j - start;

        loop {
            let rle = len.min(255);
            rle1.push(rle as u8);

            len -= rle;

            if rle != 255 {
                break;
            }
        }

        i += 1;
    }

    let mut rle2 = Vec::new();
    j = 0;
    let mut last = -1;

    while j < rle1.len() {
        let curr = rle1[j];
        j += 1;

        rle2.push(curr);

        if i32::from(curr) == last {
            let start = j;
            let mut len = 0;

            while j < rle1.len() && i32::from(rle1[j]) == last && len < 255 {
                j += 1;
                len = j - start;
            }

            rle2.push(len as u8);
        } else {
            last = i32::from(curr);
        }
    }

    writer.write_all(&rle2)?;

    Ok(())
}

fn encode_length<W>(
    writer: &mut W,
    range_coder: &mut RangeCoder,
    models: &mut Models,
    len: usize,
) -> io::Result<()>
where
    W: Write,
{
    let n = u32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    models.len[0].encode(writer, range_coder, (n & 0xff) as u8)?;
    models.len[1].encode(writer, range_coder, ((n >> 8) & 0xff) as u8)?;
    models.len[2].encode(writer, range_coder, ((n >> 16) & 0xff) as u8)?;
    models.len[3].encode(writer, range_coder, ((n >> 24) & 0xff) as u8)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codecs::fqzcomp;

    #[test]
    fn test_encode() -> io::Result<()> {
        let data = [
            vec![0, 0, 0, 1, 1, 2, 1, 1, 0, 0],
            vec![0, 1, 2, 3, 3, 3, 3, 3, 3, 3],
            vec![2, 1, 1, 0, 0],
        ];

        let records: Vec<_> = data.iter().map(|scores| (scores.len(), false)).collect();
        let src: Vec<_> = data.into_iter().flatten().collect();

        let actual = encode(&records, &src)?;

        // Verify round-trip
        let decoded = fqzcomp::decode(&actual)?;
        assert_eq!(decoded, src);

        Ok(())
    }

    #[test]
    fn test_encode_with_do_len() -> io::Result<()> {
        let data = [
            vec![0, 0, 0, 1, 1, 2, 1, 1, 0, 0],
            vec![0, 1, 2, 3, 3, 3, 3, 3, 3, 3],
            vec![2, 1, 1, 0, 0, 0, 0, 0, 1, 1],
        ];

        let records: Vec<_> = data.iter().map(|scores| (scores.len(), false)).collect();
        let src: Vec<_> = data.into_iter().flatten().collect();

        let actual = encode(&records, &src)?;

        // Verify round-trip
        let decoded = fqzcomp::decode(&actual)?;
        assert_eq!(decoded, src);

        Ok(())
    }

    #[test]
    fn test_encode_with_reverse() -> io::Result<()> {
        let data = [
            vec![10, 20, 30, 40, 50],
            vec![15, 25, 35, 45, 55],
            vec![5, 10, 15, 20, 25],
        ];

        // Second record is reverse-strand
        let records = vec![(5, false), (5, true), (5, false)];
        let src: Vec<_> = data.into_iter().flatten().collect();

        let encoded = encode(&records, &src)?;
        let decoded = fqzcomp::decode(&encoded)?;

        // After decode, reverse-strand records should have their quality scores reversed back
        // The decoder reverses them, so the decoded output should match the original
        assert_eq!(decoded, src);

        Ok(())
    }

    #[test]
    fn test_encode_decode_round_trip_varied() -> io::Result<()> {
        // Test with various quality values
        let data = [
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            vec![10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            vec![20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
        ];

        let records: Vec<_> = data.iter().map(|scores| (scores.len(), false)).collect();
        let src: Vec<_> = data.into_iter().flatten().collect();

        let encoded = encode(&records, &src)?;
        let decoded = fqzcomp::decode(&encoded)?;
        assert_eq!(decoded, src);

        Ok(())
    }

    #[test]
    fn test_encode_with_qmap_3_distinct() -> io::Result<()> {
        // 3 distinct quality values → triggers HAVE_QMAP
        let data = [
            vec![5, 10, 15, 5, 10, 15, 5, 10, 15, 5],
            vec![10, 15, 5, 10, 15, 5, 10, 15, 5, 10],
            vec![15, 5, 10, 15, 5, 10, 15, 5, 10, 15],
        ];

        let records: Vec<_> = data.iter().map(|scores| (scores.len(), false)).collect();
        let src: Vec<_> = data.into_iter().flatten().collect();

        let encoded = encode(&records, &src)?;
        let decoded = fqzcomp::decode(&encoded)?;
        assert_eq!(decoded, src);

        Ok(())
    }

    #[test]
    fn test_encode_with_qmap_single_value() -> io::Result<()> {
        // Single distinct value → triggers HAVE_QMAP (edge case)
        let data = [
            vec![7, 7, 7, 7, 7],
            vec![7, 7, 7, 7, 7],
            vec![7, 7, 7, 7, 7],
        ];

        let records: Vec<_> = data.iter().map(|scores| (scores.len(), false)).collect();
        let src: Vec<_> = data.into_iter().flatten().collect();

        let encoded = encode(&records, &src)?;
        let decoded = fqzcomp::decode(&encoded)?;
        assert_eq!(decoded, src);

        Ok(())
    }

    #[test]
    fn test_encode_with_qmap_16_values() -> io::Result<()> {
        // Exactly 16 distinct values → boundary for HAVE_QMAP
        let data: Vec<Vec<u8>> = (0..3)
            .map(|offset| (0u8..16).map(|v| v * 4 + offset).collect())
            .collect();

        let records: Vec<_> = data.iter().map(|scores| (scores.len(), false)).collect();
        let src: Vec<_> = data.into_iter().flatten().collect();

        let encoded = encode(&records, &src)?;
        let decoded = fqzcomp::decode(&encoded)?;
        assert_eq!(decoded, src);

        Ok(())
    }

    #[test]
    fn test_encode_with_no_qmap_17_values() -> io::Result<()> {
        // 17 distinct values → no HAVE_QMAP
        let data: Vec<Vec<u8>> = (0..3)
            .map(|offset| (0u8..17).map(|v| v * 3 + offset).collect())
            .collect();

        let records: Vec<_> = data.iter().map(|scores| (scores.len(), false)).collect();
        let src: Vec<_> = data.into_iter().flatten().collect();

        let encoded = encode(&records, &src)?;
        let decoded = fqzcomp::decode(&encoded)?;
        assert_eq!(decoded, src);

        Ok(())
    }

    #[test]
    fn test_build_quality_table_non_decreasing() {
        for max_q in [1, 3, 10, 15, 50, 100, 255] {
            for q_bits in [4, 5, 8, 9] {
                let table = build_quality_table(max_q, q_bits);
                assert_eq!(table.len(), 256);
                for w in table.windows(2) {
                    assert!(
                        w[0] <= w[1],
                        "not non-decreasing: max_q={max_q}, q_bits={q_bits}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_encode_with_multi_param() -> io::Result<()> {
        // 15 short (5-byte) + 15 long (20-byte) reads → triggers MULTI_PARAM
        let mut data: Vec<Vec<u8>> = Vec::new();
        for i in 0..15 {
            data.push(vec![(i % 4) as u8; 5]);
        }
        for i in 0..15 {
            data.push(vec![(i % 6) as u8; 20]);
        }

        let records: Vec<_> = data.iter().map(|scores| (scores.len(), false)).collect();
        let src: Vec<_> = data.into_iter().flatten().collect();

        let encoded = encode(&records, &src)?;
        let decoded = fqzcomp::decode(&encoded)?;
        assert_eq!(decoded, src);

        Ok(())
    }

    #[test]
    fn test_encode_uniform_lengths_no_multi_param() -> io::Result<()> {
        // Uniform lengths → no MULTI_PARAM
        let data: Vec<Vec<u8>> = (0..20).map(|i| vec![(i % 5) as u8; 10]).collect();

        let records: Vec<_> = data.iter().map(|scores| (scores.len(), false)).collect();
        let src: Vec<_> = data.into_iter().flatten().collect();

        let encoded = encode(&records, &src)?;
        let decoded = fqzcomp::decode(&encoded)?;
        assert_eq!(decoded, src);

        Ok(())
    }

    #[test]
    fn test_encode_too_few_records_no_multi_param() -> io::Result<()> {
        // Too few records → no MULTI_PARAM
        let data = vec![vec![1, 2, 3], vec![4, 5, 6, 7, 8, 9, 10, 11, 12, 13]];

        let records: Vec<_> = data.iter().map(|scores| (scores.len(), false)).collect();
        let src: Vec<_> = data.into_iter().flatten().collect();

        let encoded = encode(&records, &src)?;
        let decoded = fqzcomp::decode(&encoded)?;
        assert_eq!(decoded, src);

        Ok(())
    }

    #[test]
    fn test_encode_multi_param_with_reverse() -> io::Result<()> {
        // MULTI_PARAM + DO_REV: reverse-strand records in both groups
        let mut data: Vec<Vec<u8>> = Vec::new();
        let mut records = Vec::new();
        for i in 0..15 {
            data.push(vec![((i * 3) % 10) as u8; 5]);
            records.push((5, i % 3 == 1)); // some reverse
        }
        for i in 0..15 {
            data.push(vec![((i * 7) % 12) as u8; 20]);
            records.push((20, i % 4 == 0)); // some reverse
        }

        let src: Vec<_> = data.into_iter().flatten().collect();

        let encoded = encode(&records, &src)?;
        let decoded = fqzcomp::decode(&encoded)?;
        assert_eq!(decoded, src);

        Ok(())
    }
}
