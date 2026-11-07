//! Tie-aware ranking utilities shared across CCC modules.

/// 1-based average ranks for finite values; NaN → rank 0.0 (excluded from ties).
pub fn average_rank_f32(v: &[f32]) -> Vec<f64> {
    let n = v.len();
    let mut indexed: Vec<(usize, f32)> = v
        .iter()
        .enumerate()
        .filter(|(_, x)| x.is_finite())
        .map(|(i, &x)| (i, x))
        .collect();
    indexed.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut ranks = vec![0.0f64; n];
    let m = indexed.len();
    let mut i = 0usize;
    while i < m {
        let mut j = i + 1;
        while j < m && indexed[j].1 == indexed[i].1 {
            j += 1;
        }
        let avg = (i + j - 1) as f64 / 2.0 + 1.0;
        for k in i..j {
            ranks[indexed[k].0] = avg;
        }
        i = j;
    }
    ranks
}

/// Pearson correlation of two equal-length f64 slices.
pub fn pearson_f64(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    if x.len() != y.len() || n < 3.0 {
        return f64::NAN;
    }
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let num: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let dx = x.iter().map(|a| (a - mx).powi(2)).sum::<f64>().sqrt();
    let dy = y.iter().map(|b| (b - my).powi(2)).sum::<f64>().sqrt();
    if dx > 1e-12 && dy > 1e-12 {
        num / (dx * dy)
    } else {
        f64::NAN
    }
}

/// Spearman ρ for two equal-length f32 slices (tie-aware average ranks).
pub fn spearman_f32(x: &[f32], y: &[f32]) -> f32 {
    if x.len() != y.len() || x.len() < 3 {
        return f32::NAN;
    }
    let rx = average_rank_f32(x);
    let ry = average_rank_f32(y);
    if rx.iter().any(|r| *r <= 0.0) || ry.iter().any(|r| *r <= 0.0) {
        return f32::NAN;
    }
    pearson_f64(&rx, &ry) as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn average_rank_handles_ties() {
        let v = vec![0.0, 0.0, 0.0, 0.2, 0.8];
        let r = average_rank_f32(&v);
        assert!((r[0] - 2.0).abs() < 1e-9);
        assert!((r[1] - 2.0).abs() < 1e-9);
        assert!((r[2] - 2.0).abs() < 1e-9);
        assert!((r[3] - 4.0).abs() < 1e-9);
        assert!((r[4] - 5.0).abs() < 1e-9);
    }

    #[test]
    fn flat_values_have_zero_spearman_trend() {
        let bins = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let vals = vec![0.0, 0.0, 0.0, 0.0, 0.0];
        let rho = spearman_f32(&bins, &vals);
        assert!(
            !rho.is_finite() || rho.abs() < 1e-6,
            "flat expression should not induce monotonic trend, got rho={rho}"
        );
    }
}
