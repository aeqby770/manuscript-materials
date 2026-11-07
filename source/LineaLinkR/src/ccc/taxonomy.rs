//! Lineage entropy curves, bifurcation detection, and LR role scoring.
//!
//! Computes per-lineage entropy tensors from binned LR intermediates,
//! detects bifurcation bins, and assigns six soft role scores
//! (instructive, transition-like, reinforcing, feedback, lateral, bystander).

use crate::ccc::LrBinResult;

// ─── Per-lineage entropy tensor `H̃_l(b)` ─────────────────────────────────────

/// Per-lineage, per-bin lineage-level entropy `H_l(b)` (before smoothing)
/// **or** `H̃_l(b)` (after MA(3) smoothing) — same in-memory layout.
///
/// Layout: `data[l * n_bins + b]` row-major, so a single lineage's curve
/// is `&data[l * n_bins ..(l+1) * n_bins]`.  Values are in `[0, 1]`
/// (already divided by `log L` according to the normalization). Bins for which the
/// denominator `sum_{ct} wsums[ct,l,b]` is zero carry `f32::NAN` so the
/// downstream C2 check can count them with `.is_finite()`.
#[derive(Debug, Clone)]
pub struct LineageEntropyTensor {
    pub data: Vec<f32>,
    pub n_lineages: usize,
    pub n_bins: usize,
}

impl LineageEntropyTensor {
    /// Build an all-zero tensor of the given shape.  Used mainly by
    /// tests and by `compute_lineage_entropy_raw`'s allocation path.
    pub fn zeros(n_lineages: usize, n_bins: usize) -> Self {
        Self {
            data: vec![0.0f32; n_lineages * n_bins],
            n_lineages,
            n_bins,
        }
    }

    /// Scalar access; panics on out-of-range indexing (caller's contract
    /// check).
    #[inline]
    pub fn get(&self, l: usize, b: usize) -> f32 {
        self.data[l * self.n_bins + b]
    }

    /// Single-lineage slice `[H_l(0), H_l(1), …, H_l(n_bins-1)]`.
    #[inline]
    pub fn lineage_slice(&self, l: usize) -> &[f32] {
        &self.data[l * self.n_bins..(l + 1) * self.n_bins]
    }

    /// Mutable single-lineage slice, for in-place smoothing passes.
    #[inline]
    pub fn lineage_slice_mut(&mut self, l: usize) -> &mut [f32] {
        &mut self.data[l * self.n_bins..(l + 1) * self.n_bins]
    }
}

/// Intermediate tensors needed by `compute_lineage_entropy_raw`.
///
/// `compute_lr_scores_bins` in `lr_scoring.rs` returns these tensors in an
/// `LrBinIntermediate`, allowing per-lineage entropy to marginalise over cell
/// types without walking the cells a second time.
///
/// Layout: `[n_types × n_lineages × n_bins]` row-major, index formula
/// `ct * (n_lineages * n_bins) + l * n_bins + b`, i.e. `ct` is the
/// outermost axis (matches `wsums` / `hwsums` in `lr_scoring.rs`).
#[derive(Debug, Clone)]
pub struct LrBinIntermediate {
    pub hwsums: Vec<f64>,
    pub wsums: Vec<f64>,
    pub n_types: usize,
    pub n_lineages: usize,
    pub n_bins: usize,
}

impl LrBinIntermediate {
    /// Total number of cells = `n_types × n_lineages × n_bins`.  Named
    /// this way (not `len`) to make the layout intent explicit; used for
    /// bounds checks before marginalising.
    #[inline]
    pub fn total(&self) -> usize {
        self.n_types * self.n_lineages * self.n_bins
    }

    /// Flat index for `(ct, l, b)` under the documented layout.
    #[inline]
    pub fn idx(&self, ct: usize, l: usize, b: usize) -> usize {
        ct * (self.n_lineages * self.n_bins) + l * self.n_bins + b
    }
}

/// Compute `H_l(b)` using option A, the supported marginalisation formula.
/// Option B is not used.
///
/// ```text
/// H_l(b) = sum_{ct} hwsums[ct,l,b] / sum_{ct} wsums[ct,l,b]
/// ```
///
/// Output shape: `LineageEntropyTensor(n_lineages × n_bins)`.  Bins where
/// the denominator is zero carry `f32::NAN` so the caller's C2
/// check can detect them with `.is_finite()`.
///
/// Pure marginalisation; O(n_types × n_lineages × n_bins), no cell walk.
pub fn compute_lineage_entropy_raw(intermediate: &LrBinIntermediate) -> LineageEntropyTensor {
    let LrBinIntermediate {
        hwsums,
        wsums,
        n_types,
        n_lineages,
        n_bins,
    } = intermediate;
    let n_lineages = *n_lineages;
    let n_bins = *n_bins;
    let n_types = *n_types;

    // Degenerate shapes — return an empty-dim tensor rather than panic, so
    // an early-returned `compute_lr_scores_bins` (empty intermediate) flows
    // through cleanly.
    if n_lineages == 0 || n_bins == 0 || hwsums.is_empty() || wsums.is_empty() {
        return LineageEntropyTensor::zeros(n_lineages, n_bins);
    }

    let mut out = LineageEntropyTensor::zeros(n_lineages, n_bins);
    let slb = n_lineages * n_bins; // stride per ct (outermost axis)
    let sb = n_bins; // stride per lineage
    for l in 0..n_lineages {
        for b in 0..n_bins {
            // Marginalise sum_{ct} on both numerator and denominator per
            // the formula H_l(b) = Σ_ct hwsums / Σ_ct wsums.
            let mut num = 0.0f64;
            let mut den = 0.0f64;
            for ct in 0..n_types {
                let idx = ct * slb + l * sb + b;
                num += hwsums[idx];
                den += wsums[idx];
            }
            // Denominator 0 yields NaN so C2 can count invalid bins via
            // `.is_finite()`. The 1e-12 threshold matches `lr_scoring.rs`
            // downstream (`w > 1e-8` style); well below any real weight.
            out.data[l * n_bins + b] = if den > 1e-12 {
                (num / den) as f32
            } else {
                f32::NAN
            };
        }
    }
    out
}

/// Apply the MA(3) smoother to produce `H̃_l(b)`.
///
/// Rules:
/// - `b = 0`: `(H_l(0) + H_l(1)) / 2`  (2-bin partial window — NOT the
///   identity, NOT reflect padding).
/// - `1 <= b <= n_bins - 2`: `(H_l(b-1) + H_l(b) + H_l(b+1)) / 3`.
/// - `b = n_bins - 1`: `(H_l(n_bins-2) + H_l(n_bins-1)) / 2`.
/// - C1 bypass: if `n_bins < 6`, return `raw.clone()` unchanged.
/// - NaN handling: any window containing a NaN input propagates NaN; the
///   caller decides how to classify the affected lineage.
///
/// The smoother is non-monotone by design and allows nested bifurcation
/// patterns.
pub fn smooth_ma3(raw: &LineageEntropyTensor) -> LineageEntropyTensor {
    let n_lineages = raw.n_lineages;
    let n_bins = raw.n_bins;

    //  C1 bypass: the MA(3) smoother would eat >= 50% of the curve
    // when n_bins < 6, so the spec mandates returning the raw tensor
    // unchanged. `clone()` keeps the caller's ownership model stable.
    if n_bins < 6 {
        return raw.clone();
    }

    let mut out = LineageEntropyTensor::zeros(n_lineages, n_bins);
    for l in 0..n_lineages {
        let src_start = l * n_bins;
        // Immutable src slice; separate mutable dst slice — avoids
        // overlapping borrows the compiler would reject.
        let src_end = src_start + n_bins;
        let dst_end = src_end;

        // b = 0: 2-bin partial window (src[0] + src[1]) / 2.
        // NaN automatically propagates through `+` and `/` in IEEE-754,
        // so no explicit is_nan() check is needed.
        out.data[src_start] = (raw.data[src_start] + raw.data[src_start + 1]) * 0.5;

        // 1 ≤ b ≤ n_bins - 2: full 3-bin window.
        for b in 1..(n_bins - 1) {
            let i = src_start + b;
            out.data[i] = (raw.data[i - 1] + raw.data[i] + raw.data[i + 1]) / 3.0;
        }

        // b = n_bins - 1: 2-bin partial window (src[n-2] + src[n-1]) / 2.
        out.data[dst_end - 1] = (raw.data[src_end - 2] + raw.data[src_end - 1]) * 0.5;
    }
    out
}

// ─── Bifurcation detection and boundary tags ───────────────────────────────

/// Why a given lineage's `b_star` ended up the way it did.
///
/// The `Crossing` variant records a valid `b_star = Some(b)` from the
/// continuous-two-bin rule.
///
/// The remaining variants encode edge conditions (C2, C3, C4).
/// Only `C4BoundaryOnset` leaves `b_star = Some(b)`; C2/C3 both force
/// `b_star = None` and downgrade the lineage to "time-only" downstream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BifurcationReason {
    /// Normal case: `H̃_l(b) < h_fork` AND `H̃_l(b+1) < h_fork`.
    Crossing,

    ///  C2 — "curve not informative":
    /// `|{b : H_l(b) is finite}| < 4`, or the finite-portion
    /// `max - min < 0.05`.  `b_star = None` ⇒ lineage becomes time-only.
    C2NotInformative,

    ///  C3 — "curve informative but never enters low-entropy
    /// regime stably": `H̃_l(b)` never satisfies the continuous-2-bin
    /// rule on any interior `b`.  `b_star = None` ⇒ lineage time-only.
    C3NoLowEntropy,

    ///  C4 — boundary-onset: crossing happens at `b <= 1`.
    /// `b_star = Some(b)` but pre-fork taxonomy terms (`instructive`,
    /// `lateral`) MUST be disabled for this lineage downstream.
    /// `reinforcing` / `feedback` / `bystander` are unaffected.
    C4BoundaryOnset,
}

/// Result of `find_bifurcation_bins` search for one lineage. `b_star = None`
/// means the lineage is downgraded to time-only under C2 or C3.
/// When `reason == C4BoundaryOnset`, `b_star` is reported but downstream
/// taxonomy inference gates pre-fork terms under C4.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BifurcationResult {
    pub b_star: Option<usize>,
    pub reason: BifurcationReason,
}

impl BifurcationResult {
    /// Convenience: does this result allow pre-fork taxonomy terms
    /// (`instructive`, `lateral`) to be computed on this lineage?
    /// C4 and C2/C3 all disable them; only a clean
    /// `Crossing` allows them.
    #[inline]
    pub fn allows_prefork_terms(&self) -> bool {
        matches!(self.reason, BifurcationReason::Crossing)
    }

    /// Convenience: does this result allow post-fork taxonomy terms
    /// (`reinforcing`, `feedback`) to be computed? These
    /// are gated on `b_star.is_some()` — i.e. they work for both
    /// `Crossing` and `C4BoundaryOnset`, but not for C2/C3.
    #[inline]
    pub fn allows_postfork_terms(&self) -> bool {
        self.b_star.is_some()
    }
}

/// Find `b*(l)` for every lineage and record edge-case tags in `reason`.
///
/// Rule:
/// ```text
/// b*(l) = min { b ∈ [0, n_bins - 3] :
///                 H̃_l(b) < h_fork  ∧  H̃_l(b+1) < h_fork }
/// ```
///
/// Default `h_fork = 0.5` is a scale-natural choice. Caller is expected to
/// pass the smoothed tensor (C1-aware,
/// unsmoothed when `n_bins < 6`) — this fn does NOT re-apply MA(3).
///
/// Returns `Vec<BifurcationResult>` length `n_lineages`.
///
/// The caller applies any trajectory-level downgrade to time-only after
/// checking whether every lineage returned C2 or C3.
pub fn find_bifurcation_bins(
    h_tilde: &LineageEntropyTensor,
    h_fork: f32,
) -> Vec<BifurcationResult> {
    let n_lineages = h_tilde.n_lineages;
    let n_bins = h_tilde.n_bins;
    let mut out = Vec::with_capacity(n_lineages);

    for l in 0..n_lineages {
        let curve = h_tilde.lineage_slice(l);

        // ── C2: curve not informative ─────────────────────────────
        // Rule: n_finite < 4 OR (finite max - finite min) < 0.05.
        // Rationale: a curve with fewer than 4 finite bins, or an
        // essentially flat curve, carries no bifurcation information
        // and the caller downgrades this lineage to time-only.
        let mut n_finite: usize = 0;
        let mut fmin = f32::INFINITY;
        let mut fmax = f32::NEG_INFINITY;
        for &v in curve {
            if v.is_finite() {
                n_finite += 1;
                if v < fmin {
                    fmin = v;
                }
                if v > fmax {
                    fmax = v;
                }
            }
        }
        if n_finite < 4 || (fmax - fmin) < 0.05 {
            out.push(BifurcationResult {
                b_star: None,
                reason: BifurcationReason::C2NotInformative,
            });
            continue;
        }

        // ── Continuous two-bin crossing ───────────────────────────
        // Scan b ∈ [0, n_bins - 3] (so that b+1 ≤ n_bins - 2, matching
        // the internal-bin constraint). First `b` whose own and
        // next bin both strictly fall below `h_fork` wins.  NaN values
        // compare as false under `<`, so a NaN at b or b+1 cannot
        // trigger a false crossing.
        let mut found_b: Option<usize> = None;
        if n_bins >= 3 {
            let last_b = n_bins - 3;
            for b in 0..=last_b {
                if curve[b] < h_fork && curve[b + 1] < h_fork {
                    found_b = Some(b);
                    break;
                }
            }
        }

        // ── Classify the result by C3/C4 ───────────────────────────
        let reason = match found_b {
            // C4 boundary onset: a crossing in the first two bins leaves no
            // separable pre-fork segment. `allows_prefork_terms` disables the
            // corresponding taxonomy terms.
            Some(b) if b <= 1 => BifurcationReason::C4BoundaryOnset,
            // Clean crossing: the normal case, both pre- and post-fork
            // taxonomy terms are allowed.
            Some(_) => BifurcationReason::Crossing,
            //  C3 — informative curve but never enters low-entropy
            // regime stably: lineage gets `b_star = None` and is
            // downgraded to time-only alongside C2.
            None => BifurcationReason::C3NoLowEntropy,
        };
        out.push(BifurcationResult {
            b_star: found_b,
            reason,
        });
    }
    out
}

// ─── Role scoring tensor ────────────────────────────────────────────────────

/// Named axis for the 6-score role tensor. The order follows the three-tier
/// structure: main (instructive / transition-like / reinforcing / feedback),
/// weak (lateral), residual (bystander).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Role {
    Instructive,
    TransitionLike,
    Reinforcing,
    Feedback,
    Lateral,
    Bystander,
}

/// The six role scores.
///
/// Scores are not normalized: `sum_{role} != 1` in general.
/// - Scores are not strictly mutually exclusive; multiple non-bystander
///   terms can be non-zero on the same row.
/// - Scores do NOT claim probabilistic calibration; they are a "soft
///   role scoring system", not a classifier.
#[derive(Debug, Clone, Copy, Default)]
pub struct RoleScores {
    /// `spec(l) · 1[b < b*(l)] · σ(Δτ/τ_0)` — pre-fork, forward direction.
    /// Main tier (strong mechanism claim).
    pub instructive: f32,

    /// `spec(l) · 1[b = b*(l)]` — explicit commitment-transition bin.
    /// This prevents the inferred boundary itself from being forced into the
    /// residual bystander category.
    pub transition_like: f32,

    /// `spec(l) · 1[b > b*(l)] · σ(Δτ/τ_0)` — post-fork, forward.
    /// Main tier.
    pub reinforcing: f32,

    /// `spec(l) · 1[b > b*(l)] · σ(-Δτ/τ_0)` — post-fork, backward.
    /// Main tier.
    pub feedback: f32,

    /// `max_{l' != l} spec(l') · 1[b < b*(l)] · σ(-Δτ/τ_0)` — weak tier
    /// (for explanation only; it is not used for LR ranking).
    pub lateral: f32,

    /// `1 - max(instructive, transition_like, reinforcing, feedback, lateral)` — residual
    /// QC only; it does not participate in LR ranking or cross-pair FDR.
    pub bystander: f32,
}

impl RoleScores {
    /// Look up a single score by `Role` axis.  Useful for generic
    /// reducers (e.g. "max over main tier") in downstream consumers.
    #[inline]
    pub fn get(&self, r: Role) -> f32 {
        match r {
            Role::Instructive => self.instructive,
            Role::TransitionLike => self.transition_like,
            Role::Reinforcing => self.reinforcing,
            Role::Feedback => self.feedback,
            Role::Lateral => self.lateral,
            Role::Bystander => self.bystander,
        }
    }

    /// `max(instructive, transition_like, reinforcing, feedback, lateral)` —
    /// the value used in `bystander = 1 - this`.
    #[inline]
    pub fn max_non_bystander(&self) -> f32 {
        self.instructive
            .max(self.transition_like)
            .max(self.reinforcing)
            .max(self.feedback)
            .max(self.lateral)
    }
}

/// Parameters for the role-score computation.
#[derive(Debug, Clone, Copy)]
pub struct RoleScoreParams {
    /// `τ_0` is the sigmoid scale for the Δτ direction term. A typical value
    /// is 10% of the lineage pseudotime standard deviation. The R wrapper
    /// computes and passes this value; `tau_0 <= 0` is invalid input.
    pub tau_0: f32,

    /// `h_fork` used to derive the `BifurcationResult` consumed here.
    /// `compute_role_scores_one` does not recompute the bifurcation.
    pub h_fork: f32,

    /// **reliability gate** — minimum effective (weighted) cell count
    /// `n_eff_cells = min(source_cluster, target_cluster)` required for
    /// a row's mean-based LR score to be considered reliable enough to
    /// participate in the taxonomy. Rows with strictly less
    /// `n_eff_cells` collapse to `(0, 0, 0, 0, 0, 1)` (bystander = 1) the
    /// same way C2/C3 lineages do, but at row level rather than
    /// lineage-level granularity.
    ///
    /// Default `5.0` matches the field-level threshold documented at
    /// `src/ccc/lr_scoring.rs` (`LrBinResult::n_eff_cells`)
    /// ("Values < 5 indicate unreliable mean estimates").
    ///
    /// Set to `0.0` to disable the gate (e.g. unit tests that
    /// isolate formula behavior from row-level reliability).
    pub min_n_eff_cells: f32,

    /// Half-width (in bins) of the bifurcation transition zone.
    /// `transition_ind = 1` when `|b - b*(l)| <= transition_window`.
    /// `0` recovers exact-bin matching (`b == b*(l)` only).
    pub transition_window: u32,
}

impl Default for RoleScoreParams {
    fn default() -> Self {
        Self {
            tau_0: 1.0,
            h_fork: 0.5,
            min_n_eff_cells: 5.0,
            transition_window: 0,
        }
    }
}

/// Compute the six role scores for a single `LrBinResult`.
///
/// Arguments beyond `bin`:
///
/// - `spec_self`: `spec(p, s, t, b, l)` specificity value on this row's
///   own lineage. Caller supplies it according to the role-scoring equations.
/// - `spec_max_other`: `max_{l' != l} spec(p, s, t, b, l')` — required
///   for the `lateral` term only. The caller supplies this value; this
///   function does not reduce across rows.
/// - `delta_tau`: Δτ = target_pt - source_pt per row (positive when the
///   receiver pseudotime exceeds the sender, i.e. sender leads).
/// - `b_star`:  result for this row's lineage, including the
///   C2/C3/C4 reason tag.
/// - `params`: `τ_0` and the recorded `h_fork` setting; the row-level fork
///   position is supplied through `b_star`.
///
/// C2/C3 (`b_star.b_star == None`) return `(0,0,0,0,0,1)`. Under C4,
/// `instructive` and `lateral` are zero even when their indicators would fire.
pub fn compute_role_scores_one(
    bin: &LrBinResult,
    spec_self: f32,
    spec_max_other: f32,
    delta_tau: f32,
    b_star: BifurcationResult,
    params: &RoleScoreParams,
) -> RoleScores {
    // ── C2/C3 early exit ─────────────────────────────────────────────
    // `b_star.b_star = None` ⇒ this lineage has been downgraded to
    // time-only; every row on this lineage receives `(0,0,0,0,0,1)`.
    // The natural fall-through produces the same values. The explicit return
    // keeps the C2/C3 semantics visible at the top of the function.
    if b_star.b_star.is_none() {
        return RoleScores {
            instructive: 0.0,
            transition_like: 0.0,
            reinforcing: 0.0,
            feedback: 0.0,
            lateral: 0.0,
            bystander: 1.0,
        };
    }

    // ── row-level reliability gate ──────────────────────────────────────
    // If `n_eff_cells < min_n_eff_cells`, the mean-based LR score is treated
    // as unreliable and reduced to bystander, matching the C2/C3 lineage gate.
    // Non-finite `n_eff_cells` values are not below the threshold; callers
    // should provide finite values. `min_n_eff_cells = 0.0` disables this gate.
    if bin.n_eff_cells < params.min_n_eff_cells {
        return RoleScores {
            instructive: 0.0,
            transition_like: 0.0,
            reinforcing: 0.0,
            feedback: 0.0,
            lateral: 0.0,
            bystander: 1.0,
        };
    }

    let bs = b_star.b_star.unwrap();
    let b = bin.bin;

    // ── Δτ-direction sigmoid ──────────────────────────────────────────
    // τ_0 <= 0 is a caller error. The R wrapper validates this input before
    // calling Rust.
    debug_assert!(params.tau_0 > 0.0, "RoleScoreParams.tau_0 must be > 0");
    let dt_scaled = delta_tau / params.tau_0;
    let sig_fwd = sigmoid(dt_scaled); // σ( Δτ / τ_0)
    let sig_bwd = sigmoid(-dt_scaled); // σ(-Δτ / τ_0)

    // ── Pre/post indicators with C4 gating ────────────────────────────
    //
    // Pre-fork indicator `1[b < b*(l)]` fires ONLY for the clean
    // `Crossing` reason (C4 explicitly disables pre-fork terms
    // on this lineage).  `BifurcationResult::allows_prefork_terms`
    // encodes this directly: it matches on `Crossing` alone.
    //
    // Post-fork indicator `1[b > b*(l)]` is allowed whenever
    // `b_star.is_some()`, i.e. both `Crossing` and `C4BoundaryOnset`
    // Reinforcing, feedback, and bystander scores are unaffected by this gate.
    let pre_ind: f32 = if (b < bs) && b_star.allows_prefork_terms() {
        1.0
    } else {
        0.0
    };
    let post_ind: f32 = if (b > bs) && b_star.allows_postfork_terms() {
        1.0
    } else {
        0.0
    };
    let dist_bins = if b >= bs { b - bs } else { bs - b };
    let tw = params.transition_window as usize;
    let transition_ind: f32 = if dist_bins <= tw { 1.0 } else { 0.0 };

    // ── six formulas, including the explicit transition bin ─────────
    //   S_instructive = spec(l) · 1[b < b*(l)] · σ( Δτ / τ_0)  (pre, fwd)
    //   S_reinforcing = spec(l) · 1[b > b*(l)] · σ( Δτ / τ_0)  (post, fwd)
    //   S_feedback    = spec(l) · 1[b > b*(l)] · σ(-Δτ / τ_0)  (post, bwd)
    //   S_lateral     = max_{l'≠l} spec(l') · 1[b < b*(l)] · σ(-Δτ / τ_0)
    //                                                          (pre, bwd, opposing)
    //   S_transition  = spec(l) · 1[|b - b*(l)| <= transition_window]
    //   S_bystander   = 1 - max(instructive, transition, reinforcing, feedback, lateral)
    let instructive = spec_self * pre_ind * sig_fwd;
    let transition_like = spec_self * transition_ind;
    let reinforcing = spec_self * post_ind * sig_fwd;
    let feedback = spec_self * post_ind * sig_bwd;
    let lateral = spec_max_other * pre_ind * sig_bwd;

    let max_main = instructive
        .max(transition_like)
        .max(reinforcing)
        .max(feedback)
        .max(lateral);
    let bystander = if [instructive, transition_like, reinforcing, feedback, lateral]
        .iter()
        .any(|score| score.is_nan())
    {
        f32::NAN
    } else {
        1.0 - max_main
    };

    RoleScores {
        instructive,
        transition_like,
        reinforcing,
        feedback,
        lateral,
        bystander,
    }
}

/// Standard logistic sigmoid.  Private taxonomy helper — the `-x.exp()`
/// form is stable for large positive `x` (no overflow before negation),
/// and `NaN` input propagates to `NaN` output as expected.
#[inline]
fn sigmoid(x: f32) -> f32 {
    1.0 / (1.0 + (-x).exp())
}

#[cfg(test)]
mod tests {
    //! Unit tests for the role-score data structures and implemented taxonomy
    //! kernels.

    use super::*;

    #[test]
    fn lineage_entropy_tensor_layout_round_trip() {
        let mut t = LineageEntropyTensor::zeros(3, 4);
        // Distinct values per (l, b) test the documented
        // `l * n_bins + b` layout.
        for l in 0..3 {
            for b in 0..4 {
                t.data[l * 4 + b] = (10 * l + b) as f32;
            }
        }
        assert_eq!(t.get(0, 0), 0.0);
        assert_eq!(t.get(1, 2), 12.0);
        assert_eq!(t.get(2, 3), 23.0);
        assert_eq!(t.lineage_slice(1), &[10.0, 11.0, 12.0, 13.0]);
        t.lineage_slice_mut(0)[2] = 99.0;
        assert_eq!(t.get(0, 2), 99.0);
    }

    #[test]
    fn lr_bin_intermediate_idx_matches_documented_layout() {
        // 2 types × 3 lineages × 4 bins = 24 cells.
        let mut inter = LrBinIntermediate {
            hwsums: vec![0.0; 24],
            wsums: vec![0.0; 24],
            n_types: 2,
            n_lineages: 3,
            n_bins: 4,
        };
        assert_eq!(inter.total(), 24);
        // ct outermost, then l, then b.
        assert_eq!(inter.idx(0, 0, 0), 0);
        assert_eq!(inter.idx(0, 0, 3), 3);
        assert_eq!(inter.idx(0, 1, 0), 4);
        assert_eq!(inter.idx(0, 2, 3), 11);
        assert_eq!(inter.idx(1, 0, 0), 12);
        assert_eq!(inter.idx(1, 2, 3), 23);
        // Write-read round trip to confirm no off-by-one in idx().
        let i = inter.idx(1, 1, 2);
        inter.hwsums[i] = 7.5;
        assert_eq!(inter.hwsums[inter.idx(1, 1, 2)], 7.5);
    }

    #[test]
    fn bifurcation_result_gating_helpers_match_prose() {
        //  C4 table: post-fork terms allowed, pre-fork disabled.
        let r_c4 = BifurcationResult {
            b_star: Some(1),
            reason: BifurcationReason::C4BoundaryOnset,
        };
        assert!(!r_c4.allows_prefork_terms());
        assert!(r_c4.allows_postfork_terms());

        // Normal crossing: both allowed.
        let r_ok = BifurcationResult {
            b_star: Some(5),
            reason: BifurcationReason::Crossing,
        };
        assert!(r_ok.allows_prefork_terms());
        assert!(r_ok.allows_postfork_terms());

        // C2 / C3: neither allowed (lineage is time-only).
        for reason in [
            BifurcationReason::C2NotInformative,
            BifurcationReason::C3NoLowEntropy,
        ] {
            let r = BifurcationResult {
                b_star: None,
                reason,
            };
            assert!(!r.allows_prefork_terms());
            assert!(!r.allows_postfork_terms());
        }
    }

    #[test]
    fn role_scores_get_and_max_non_bystander() {
        let s = RoleScores {
            instructive: 0.4,
            transition_like: 0.2,
            reinforcing: 0.7,
            feedback: 0.1,
            lateral: 0.3,
            bystander: 0.3, // = 1 - max(0.4, 0.7, 0.1, 0.3)
        };
        assert_eq!(s.get(Role::Instructive), 0.4);
        assert_eq!(s.get(Role::TransitionLike), 0.2);
        assert_eq!(s.get(Role::Reinforcing), 0.7);
        assert_eq!(s.get(Role::Feedback), 0.1);
        assert_eq!(s.get(Role::Lateral), 0.3);
        assert_eq!(s.get(Role::Bystander), 0.3);
        assert!((s.max_non_bystander() - 0.7).abs() < 1e-6);
    }

    // ─── Per-lineage entropy tests ──────────────────────────────────────────

    /// Basic marginalisation: H_l(b) = Σ_ct hwsums[ct,l,b] / Σ_ct wsums[ct,l,b].
    /// Hand-constructed 2-type × 2-lineage × 3-bin fixture for direct
    /// comparison with paper-and-pencil values.
    #[test]
    fn compute_lineage_entropy_raw_basic_marginalisation() {
        // Layout convention: ct * (n_lineages * n_bins) + l * n_bins + b.
        //   n_types = 2, n_lineages = 2, n_bins = 3, total 12 cells.
        //
        // Design the numbers so that (l=0, b=0) has ct-sum 4/10 = 0.4,
        //   (l=0, b=1) has 6/12 = 0.5, (l=1, b=2) has 9/18 = 0.5, etc.
        let hwsums: Vec<f64> = vec![
            // ct = 0
            1.0, 2.0, 3.0, // l=0  (b=0,1,2)
            4.0, 5.0, 6.0, // l=1
            // ct = 1
            3.0, 4.0, 5.0, // l=0
            2.0, 3.0, 4.0, // l=1
        ];
        let wsums: Vec<f64> = vec![
            // ct = 0
            4.0, 6.0, 10.0, 5.0, 5.0, 5.0, // ct = 1
            6.0, 6.0, 10.0, 5.0, 7.0, 13.0,
        ];
        let inter = LrBinIntermediate {
            hwsums,
            wsums,
            n_types: 2,
            n_lineages: 2,
            n_bins: 3,
        };
        let out = compute_lineage_entropy_raw(&inter);
        assert_eq!(out.n_lineages, 2);
        assert_eq!(out.n_bins, 3);

        // l=0, b=0:  num = 1+3 = 4,  den = 4+6 = 10   → 0.4
        // l=0, b=1:  num = 2+4 = 6,  den = 6+6 = 12   → 0.5
        // l=0, b=2:  num = 3+5 = 8,  den = 10+10=20  → 0.4
        // l=1, b=0:  num = 4+2 = 6,  den = 5+5 = 10   → 0.6
        // l=1, b=1:  num = 5+3 = 8,  den = 5+7 = 12   → 0.6666...
        // l=1, b=2:  num = 6+4 = 10, den = 5+13=18   → 0.5555...
        let expected: [f32; 6] = [0.4, 0.5, 0.4, 0.6, 8.0 / 12.0, 10.0 / 18.0];
        for (i, &e) in expected.iter().enumerate() {
            assert!(
                (out.data[i] - e).abs() < 1e-5,
                "idx {i}: expected {e}, got {}",
                out.data[i]
            );
        }
    }

    ///  C2 invariance: a (l, b) cell whose denominator Σ_ct wsums
    /// is zero MUST come back as NaN, so the caller can count invalid
    /// bins with `.is_finite()` instead of probing numerator directly.
    #[test]
    fn compute_lineage_entropy_raw_nan_on_zero_denominator() {
        // 1 type × 1 lineage × 2 bins.  bin 0 has wsums = 0 → NaN expected.
        let inter = LrBinIntermediate {
            hwsums: vec![0.0, 1.0],
            wsums: vec![0.0, 2.0], // bin 0: zero denom
            n_types: 1,
            n_lineages: 1,
            n_bins: 2,
        };
        let out = compute_lineage_entropy_raw(&inter);
        assert!(
            out.data[0].is_nan(),
            "zero-denominator bin must be NaN, got {}",
            out.data[0]
        );
        assert!((out.data[1] - 0.5).abs() < 1e-6);
    }

    ///  C1 bypass: MA(3) on a curve shorter than 6 bins returns the
    /// raw tensor unchanged — any partial window would eat ≥ 50 % of the
    /// curve, so the smoother is skipped entirely.
    #[test]
    fn smooth_ma3_bypass_when_n_bins_below_6() {
        // n_bins = 5 → bypass.
        let raw = LineageEntropyTensor {
            data: vec![0.1, 0.2, 0.3, 0.4, 0.5],
            n_lineages: 1,
            n_bins: 5,
        };
        let smoothed = smooth_ma3(&raw);
        assert_eq!(
            smoothed.data, raw.data,
            "n_bins < 6 must return the raw tensor bit-identical"
        );
    }

    ///  boundary rule: at b = 0 and b = n_bins - 1 the window is a
    /// 2-bin partial mean, NOT a 3-bin window with identity padding and
    /// NOT a reflect / constant extrapolation.
    #[test]
    fn smooth_ma3_2bin_partial_windows_at_boundaries() {
        // n_bins = 6, single lineage, monotone increasing values so the
        // math is trivial to verify.  H = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0].
        let raw = LineageEntropyTensor {
            data: vec![0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            n_lineages: 1,
            n_bins: 6,
        };
        let out = smooth_ma3(&raw);
        // b = 0: (0.0 + 0.2) / 2 = 0.1
        assert!(
            (out.data[0] - 0.1).abs() < 1e-6,
            "boundary b=0 expected 0.1, got {}",
            out.data[0]
        );
        // b = 5 (= n_bins - 1): (0.8 + 1.0) / 2 = 0.9
        assert!(
            (out.data[5] - 0.9).abs() < 1e-6,
            "boundary b=n-1 expected 0.9, got {}",
            out.data[5]
        );
    }

    ///  interior rule: 3-bin window `(H[b-1] + H[b] + H[b+1]) / 3`
    /// for every `1 <= b <= n_bins - 2`. The monotone fixture permits direct
    /// checking of each interior bin.
    #[test]
    fn smooth_ma3_3bin_window_in_interior() {
        let raw = LineageEntropyTensor {
            data: vec![0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
            n_lineages: 1,
            n_bins: 6,
        };
        let out = smooth_ma3(&raw);
        // b = 1:  (0.0 + 0.2 + 0.4) / 3 = 0.2
        // b = 2:  (0.2 + 0.4 + 0.6) / 3 = 0.4
        // b = 3:  (0.4 + 0.6 + 0.8) / 3 = 0.6
        // b = 4:  (0.6 + 0.8 + 1.0) / 3 = 0.8
        for &(b, e) in &[(1usize, 0.2f32), (2, 0.4), (3, 0.6), (4, 0.8)] {
            assert!(
                (out.data[b] - e).abs() < 1e-6,
                "interior b={b} expected {e}, got {}",
                out.data[b]
            );
        }
    }

    /// NaN propagation: any 3-bin or 2-bin window that contains a NaN input
    /// must emit NaN. `+` and `/` propagate NaN in IEEE-754; the implementation
    /// must not short-circuit this propagation.
    #[test]
    fn smooth_ma3_nan_propagates_within_window() {
        // n_bins = 6, put NaN at b = 2.  Then:
        //   b = 0  window [0,1]      → finite  (no NaN)
        //   b = 1  window [0,1,2]    → NaN     (contains b=2)
        //   b = 2  window [1,2,3]    → NaN
        //   b = 3  window [2,3,4]    → NaN
        //   b = 4  window [3,4,5]    → finite
        //   b = 5  window [4,5]      → finite
        let raw = LineageEntropyTensor {
            data: vec![0.0, 0.2, f32::NAN, 0.6, 0.8, 1.0],
            n_lineages: 1,
            n_bins: 6,
        };
        let out = smooth_ma3(&raw);
        assert!(out.data[0].is_finite());
        assert!(out.data[1].is_nan());
        assert!(out.data[2].is_nan());
        assert!(out.data[3].is_nan());
        assert!(out.data[4].is_finite());
        assert!(out.data[5].is_finite());
    }

    // ─── Bifurcation detection tests (C2/C3/C4) ────────────────────────────

    /// Helper: single-lineage tensor from a per-bin slice, so the bifurcation detection
    /// tests below read like specifications.
    fn single_lineage_tensor(curve: Vec<f32>) -> LineageEntropyTensor {
        let n_bins = curve.len();
        LineageEntropyTensor {
            data: curve,
            n_lineages: 1,
            n_bins,
        }
    }

    /// A smooth monotone-decreasing curve crosses
    /// `h_fork = 0.5` at an interior bin, `b_star` is the first such
    /// `b`, `reason = Crossing`, pre- and post-fork terms both allowed.
    #[test]
    fn find_bifurcation_bins_normal_crossing() {
        // 8 bins, H̃ = [0.9, 0.85, 0.7, 0.55, 0.45, 0.40, 0.35, 0.30].
        // Strict < 0.5: b=4 (0.45) and b=5 (0.40) both pass.
        // b=3 (0.55) fails strict <; so first crossing at b=4.
        let tensor = single_lineage_tensor(vec![0.9, 0.85, 0.7, 0.55, 0.45, 0.40, 0.35, 0.30]);
        let results = find_bifurcation_bins(&tensor, 0.5);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].b_star, Some(4));
        assert_eq!(results[0].reason, BifurcationReason::Crossing);
        assert!(results[0].allows_prefork_terms());
        assert!(results[0].allows_postfork_terms());
    }

    ///  C2 via n_finite < 4: only 3 of 8 bins are finite, so the
    /// lineage is flagged `c2_not_informative` regardless of whether
    /// those 3 values look like a crossing.
    #[test]
    fn find_bifurcation_bins_c2_too_few_finite() {
        let nan = f32::NAN;
        let tensor = single_lineage_tensor(vec![0.3, 0.2, 0.1, nan, nan, nan, nan, nan]);
        let results = find_bifurcation_bins(&tensor, 0.5);
        assert_eq!(results[0].b_star, None);
        assert_eq!(results[0].reason, BifurcationReason::C2NotInformative);
        assert!(!results[0].allows_prefork_terms());
        assert!(!results[0].allows_postfork_terms());
    }

    ///  C2 via flat curve: 8 finite bins but finite max - min
    /// < 0.05 (strict less-than), so the lineage is flagged
    /// `c2_not_informative` even though some bins sit below `h_fork`.
    #[test]
    fn find_bifurcation_bins_c2_flat_curve() {
        // Spread = 0.49 - 0.46 = 0.03, strictly < 0.05 threshold.
        let tensor = single_lineage_tensor(vec![0.46, 0.47, 0.48, 0.49, 0.49, 0.48, 0.47, 0.46]);
        let results = find_bifurcation_bins(&tensor, 0.5);
        assert_eq!(results[0].b_star, None);
        assert_eq!(results[0].reason, BifurcationReason::C2NotInformative);
    }

    ///  C3: curve is informative (8 finite, spread 0.6) but never
    /// dips two consecutive bins below `h_fork = 0.5`.
    #[test]
    fn find_bifurcation_bins_c3_never_crosses() {
        // Every bin stays >= 0.5 (so strict < fails everywhere).
        let tensor = single_lineage_tensor(vec![0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60]);
        let results = find_bifurcation_bins(&tensor, 0.5);
        assert_eq!(results[0].b_star, None);
        assert_eq!(results[0].reason, BifurcationReason::C3NoLowEntropy);
        assert!(!results[0].allows_prefork_terms());
        assert!(!results[0].allows_postfork_terms());
    }

    ///  C4: crossing happens at b = 0 or b = 1.  `b_star` is
    /// reported, but `allows_prefork_terms()` must return false so
    /// downstream taxonomy gates the `instructive` / `lateral` terms.
    /// `allows_postfork_terms()` stays true — `reinforcing` / `feedback`
    /// are still well-defined.
    #[test]
    fn find_bifurcation_bins_c4_boundary_onset() {
        // Crossing at b = 0 (both 0.3 and 0.3 < 0.5).  Informative curve
        // (spread >= 0.05, > 4 finite).
        let tensor = single_lineage_tensor(vec![0.30, 0.30, 0.40, 0.55, 0.60, 0.65, 0.70, 0.80]);
        let results = find_bifurcation_bins(&tensor, 0.5);
        assert_eq!(results[0].b_star, Some(0));
        assert_eq!(results[0].reason, BifurcationReason::C4BoundaryOnset);
        assert!(
            !results[0].allows_prefork_terms(),
            "C4 must disable pre-fork taxonomy terms"
        );
        assert!(
            results[0].allows_postfork_terms(),
            "C4 still allows post-fork taxonomy terms"
        );

        // Crossing at b = 1 also counts as C4.
        let tensor = single_lineage_tensor(vec![0.80, 0.30, 0.30, 0.55, 0.60, 0.65, 0.70, 0.80]);
        let results = find_bifurcation_bins(&tensor, 0.5);
        assert_eq!(results[0].b_star, Some(1));
        assert_eq!(results[0].reason, BifurcationReason::C4BoundaryOnset);
    }

    ///  continuous-2-bin rule rejects a single-bin dip.  A curve
    /// that briefly drops below `h_fork` for one bin but comes back up
    /// does NOT constitute a crossing — the bifurcation is not stable.
    #[test]
    fn find_bifurcation_bins_rejects_single_bin_dip() {
        // Only b=2 is below 0.5; b=3 bounces back to 0.6.
        // No two consecutive bins < 0.5 anywhere on the curve, so C3.
        let tensor = single_lineage_tensor(vec![0.9, 0.8, 0.4, 0.6, 0.55, 0.52, 0.51, 0.55]);
        let results = find_bifurcation_bins(&tensor, 0.5);
        assert_eq!(results[0].b_star, None);
        assert_eq!(results[0].reason, BifurcationReason::C3NoLowEntropy);
    }

    /// NaN at a candidate bin must not trigger a false crossing because
    /// `NaN < h_fork` is false under IEEE-754. An MA(3) window containing NaN
    /// propagates it to that bin, which remains ineligible for crossing.
    #[test]
    fn find_bifurcation_bins_nan_does_not_trigger() {
        let nan = f32::NAN;
        // 8 bins: b=3 is NaN.  Ignoring b=3, the "two consecutive below"
        // rule first fires at b=5 (0.45 and 0.40).  Must NOT report
        // b_star = 3 just because the NaN window is "ambiguous".
        //
        // Finite count = 7 (>=4), finite spread = 0.9 - 0.30 = 0.60 (>=0.05),
        // so this is NOT C2 either.
        let tensor = single_lineage_tensor(vec![0.9, 0.85, 0.7, nan, 0.55, 0.45, 0.40, 0.30]);
        let results = find_bifurcation_bins(&tensor, 0.5);
        // b=5 wins: curve[5]=0.45 < 0.5 AND curve[6]=0.40 < 0.5.
        // b=3 (NaN) cannot win because NaN < 0.5 is false.
        assert_eq!(results[0].b_star, Some(5));
        assert_eq!(results[0].reason, BifurcationReason::Crossing);
    }

    /// Multi-lineage driver: each lineage is judged independently, and
    /// the returned Vec preserves input order.  Covers all 4 reason
    /// variants in one call to catch any accidental state leakage
    /// across lineages.
    #[test]
    fn find_bifurcation_bins_multi_lineage_independent() {
        let nan = f32::NAN;
        // 4 lineages × 8 bins, one per reason variant.
        let mut data = Vec::with_capacity(4 * 8);
        //  l=0 Crossing:   first 2-consec below 0.5 at b=4.
        data.extend_from_slice(&[0.9, 0.85, 0.7, 0.55, 0.45, 0.40, 0.35, 0.30]);
        //  l=1 C2 (flat, spread = 0.03 < 0.05 strict threshold):
        data.extend_from_slice(&[0.46, 0.47, 0.48, 0.49, 0.49, 0.48, 0.47, 0.46]);
        //  l=2 C3 (never below 0.5 two bins in a row):
        data.extend_from_slice(&[0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60]);
        //  l=3 C4 (crossing at b=0):
        data.extend_from_slice(&[0.30, 0.30, 0.40, 0.55, 0.60, 0.65, 0.70, 0.80]);
        let _ = nan; // silence unused_imports-style warning on the helper

        let tensor = LineageEntropyTensor {
            data,
            n_lineages: 4,
            n_bins: 8,
        };
        let results = find_bifurcation_bins(&tensor, 0.5);
        assert_eq!(results.len(), 4);
        assert_eq!(results[0].reason, BifurcationReason::Crossing);
        assert_eq!(results[0].b_star, Some(4));
        assert_eq!(results[1].reason, BifurcationReason::C2NotInformative);
        assert_eq!(results[1].b_star, None);
        assert_eq!(results[2].reason, BifurcationReason::C3NoLowEntropy);
        assert_eq!(results[2].b_star, None);
        assert_eq!(results[3].reason, BifurcationReason::C4BoundaryOnset);
        assert_eq!(results[3].b_star, Some(0));
    }

    /// Multi-lineage layout check: rows of the LineageEntropyTensor are
    /// smoothed independently — no cross-lineage contamination in the
    /// 3-bin window — and the per-lineage `n_bins` stride is respected.
    #[test]
    fn smooth_ma3_lineages_are_independent() {
        // 2 lineages × 6 bins.  Lineage 0 constant 0.5, lineage 1
        // constant 0.9.  After MA(3), both curves should be exactly
        // their respective constants (window average of a constant is
        // the constant).  Any cross-lineage leakage would show up as
        // a weighted average between 0.5 and 0.9.
        let raw = LineageEntropyTensor {
            data: vec![0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9],
            n_lineages: 2,
            n_bins: 6,
        };
        let out = smooth_ma3(&raw);
        for b in 0..6 {
            assert!((out.data[b] - 0.5).abs() < 1e-6, "lineage 0 bin {b}");
            assert!((out.data[6 + b] - 0.9).abs() < 1e-6, "lineage 1 bin {b}");
        }
    }

    // ─── Role scoring tests (six formulas plus C2/C3/C4 gating) ────────────

    /// Helper: construct an `LrBinResult` stub that only the `.bin` and
    /// `.lineage` fields matter for `compute_role_scores_one` (per
    /// the function reads `bin.bin` and nothing else from the struct).
    fn bin_stub(bin_idx: usize, lineage_idx: usize) -> LrBinResult {
        LrBinResult {
            pair_idx: 0,
            source: 0,
            target: 0,
            lineage: lineage_idx,
            bin: bin_idx,
            score: 0.0,
            n_eff_cells: 0.0,
            weight_mass_cells: 0.0,
            lineage_entropy: 0.0,
        }
    }

    ///  C2/C3 early exit: a `BifurcationResult` with `b_star = None`
    /// (either C2 or C3) must collapse the 6-score tuple to (0,0,0,0,0,1)
    /// regardless of spec / Δτ — the whole lineage is downgraded to
    /// time-only so no taxonomy claim is made on this row.
    #[test]
    fn compute_role_scores_one_c2_c3_early_return() {
        let bin = bin_stub(5, 0);
        // Disable the reliability gate so this test isolates the role formula.
        let params = RoleScoreParams {
            tau_0: 1.0,
            h_fork: 0.5,
            min_n_eff_cells: 0.0,
            transition_window: 0,
        };

        for reason in [
            BifurcationReason::C2NotInformative,
            BifurcationReason::C3NoLowEntropy,
        ] {
            let bif = BifurcationResult {
                b_star: None,
                reason,
            };
            // Intentionally make every other input "hot" so only the
            // early-exit path can produce (0,0,0,0,0,1).
            let s = compute_role_scores_one(
                &bin, /*spec_self*/ 0.9, /*spec_max_other*/ 0.8, /*delta_tau*/ 10.0,
                bif, &params,
            );
            assert_eq!(s.instructive, 0.0, "reason {:?}: instructive", reason);
            assert_eq!(s.reinforcing, 0.0, "reason {:?}: reinforcing", reason);
            assert_eq!(s.feedback, 0.0, "reason {:?}: feedback", reason);
            assert_eq!(s.lateral, 0.0, "reason {:?}: lateral", reason);
            assert_eq!(s.bystander, 1.0, "reason {:?}: bystander", reason);
        }
    }

    #[test]
    fn compute_role_scores_sender_leads_receiver() {
        let bin = bin_stub(2, 0);
        let bif = BifurcationResult {
            b_star: Some(5),
            reason: BifurcationReason::Crossing,
        };
        let params = RoleScoreParams {
            tau_0: 0.1,
            h_fork: 0.5,
            min_n_eff_cells: 0.0,
            transition_window: 0,
        };
        // receiver_pt - sender_pt > 0 ⇒ sender leads ⇒ higher instructive
        let sender_leads = compute_role_scores_one(&bin, 0.8, 0.0, 0.2, bif.clone(), &params);
        let receiver_leads = compute_role_scores_one(&bin, 0.8, 0.0, -0.2, bif, &params);
        assert!(
            sender_leads.instructive > receiver_leads.instructive,
            "sender_leads instructive={}, receiver_leads instructive={}",
            sender_leads.instructive,
            receiver_leads.instructive
        );
    }

    /// Clean pre-fork, forward-direction row: b < b*(l), Δτ > 0, spec
    /// non-zero.  Only `instructive` should fire on the main tier;
    /// `reinforcing` / `feedback` are zero because post-fork indicator
    /// is off; `lateral` is non-zero only if spec_max_other > 0 AND
    /// σ(-Δτ/τ_0) is non-negligible.  With τ_0 = 1.0 and Δτ = 5.0,
    /// σ(-5) ≈ 0.00669, so lateral ≈ spec_max_other * 0.00669 — small
    /// but non-zero.
    #[test]
    fn compute_role_scores_one_prefork_forward_activates_instructive() {
        let bin = bin_stub(/*b=*/ 2, 0);
        let bif = BifurcationResult {
            b_star: Some(5), // b < b*
            reason: BifurcationReason::Crossing,
        };
        let params = RoleScoreParams {
            tau_0: 1.0,
            h_fork: 0.5,
            min_n_eff_cells: 0.0,
            transition_window: 0,
        };
        let s = compute_role_scores_one(
            &bin, /*spec_self*/ 0.8,
            /*spec_max_other*/ 0.0, // zero to isolate instructive
            /*delta_tau*/ 5.0, // sigmoid fwd ≈ 0.9933
            bif, &params,
        );

        // instructive = 0.8 · 1 · σ(5) ≈ 0.8 · 0.9933 ≈ 0.7946
        assert!(
            (s.instructive - 0.8 * 0.9933).abs() < 1e-3,
            "instructive = {}",
            s.instructive
        );
        assert_eq!(s.reinforcing, 0.0);
        assert_eq!(s.feedback, 0.0);
        assert_eq!(s.lateral, 0.0); // spec_max_other = 0
                                    // bystander = 1 - max = 1 - instructive
        assert!((s.bystander - (1.0 - s.instructive)).abs() < 1e-6);
    }

    /// Clean post-fork row with Δτ = 0: σ(0) = σ(-0) = 0.5 exactly, so
    /// `reinforcing` and `feedback` must both equal spec_self * 0.5.
    /// `instructive` / `lateral` zero because b > b*.  This nails the
    /// symmetry of the  forward/backward sigmoid pair at Δτ = 0.
    #[test]
    fn compute_role_scores_one_postfork_at_zero_delta_is_symmetric() {
        let bin = bin_stub(/*b=*/ 7, 0);
        let bif = BifurcationResult {
            b_star: Some(3), // b > b*
            reason: BifurcationReason::Crossing,
        };
        let params = RoleScoreParams {
            tau_0: 1.0,
            h_fork: 0.5,
            min_n_eff_cells: 0.0,
            transition_window: 0,
        };
        let s = compute_role_scores_one(
            &bin, /*spec_self*/ 0.6, /*spec_max_other*/ 0.5, /*delta_tau*/ 0.0, bif,
            &params,
        );

        assert_eq!(s.instructive, 0.0);
        assert!((s.reinforcing - 0.3).abs() < 1e-6); // 0.6 · 1 · 0.5
        assert!((s.feedback - 0.3).abs() < 1e-6); // 0.6 · 1 · 0.5
        assert_eq!(s.lateral, 0.0); // pre indicator = 0
        assert!((s.bystander - (1.0 - 0.3)).abs() < 1e-6);
    }

    ///  C4 gating: even when `b < b*(l)` would fire the pre-fork
    /// indicator, the C4BoundaryOnset reason MUST force `instructive`
    /// and `lateral` to zero.  `reinforcing` and `feedback` (post-fork)
    /// stay enabled in C4, so they react normally on rows with b > b*.
    #[test]
    fn compute_role_scores_one_c4_gates_prefork_terms() {
        let params = RoleScoreParams {
            tau_0: 1.0,
            h_fork: 0.5,
            min_n_eff_cells: 0.0,
            transition_window: 0,
        };

        // Case A: C4 lineage, row sits in pre-fork (b < b*).  Pre terms
        // MUST be zero; post terms are also zero because b < b*.
        let bif_c4 = BifurcationResult {
            b_star: Some(1), // C4 boundary: b* ≤ 1
            reason: BifurcationReason::C4BoundaryOnset,
        };
        let bin_pre = bin_stub(/*b=*/ 0, 0); // b < b*
        let s_pre = compute_role_scores_one(
            &bin_pre, /*spec_self*/ 0.9, /*spec_max_other*/ 0.9, /*delta_tau*/ 5.0,
            bif_c4, &params,
        );
        assert_eq!(s_pre.instructive, 0.0, "C4 must zero instructive");
        assert_eq!(s_pre.lateral, 0.0, "C4 must zero lateral");
        assert_eq!(s_pre.reinforcing, 0.0);
        assert_eq!(s_pre.feedback, 0.0);
        assert_eq!(s_pre.bystander, 1.0);

        // Case B: C4 lineage with a post-fork row (b > b*). Post-fork
        // terms remain enabled.
        let bin_post = bin_stub(/*b=*/ 5, 0); // b > b*=1
        let s_post = compute_role_scores_one(
            &bin_post, /*spec_self*/ 0.9, /*spec_max_other*/ 0.9, /*delta_tau*/ 0.0,
            bif_c4, &params,
        );
        assert_eq!(s_post.instructive, 0.0);
        assert_eq!(s_post.lateral, 0.0);
        assert!((s_post.reinforcing - 0.45).abs() < 1e-6); // 0.9·1·0.5
        assert!((s_post.feedback - 0.45).abs() < 1e-6);
    }

    /// Row exactly at `b = b*(l)`: both `b < b*` and `b > b*` indicators
    /// are false; the explicit transition score claims the bifurcation bin.
    #[test]
    fn compute_role_scores_one_at_bifurcation_bin_is_transition_like() {
        let bin = bin_stub(/*b=*/ 3, 0);
        let bif = BifurcationResult {
            b_star: Some(3), // b == b*
            reason: BifurcationReason::Crossing,
        };
        let params = RoleScoreParams {
            tau_0: 1.0,
            h_fork: 0.5,
            min_n_eff_cells: 0.0,
            transition_window: 0,
        };
        let s = compute_role_scores_one(
            &bin, /*spec_self*/ 0.9, /*spec_max_other*/ 0.9, /*delta_tau*/ 3.0, bif,
            &params,
        );
        assert_eq!(s.instructive, 0.0);
        assert!((s.transition_like - 0.9).abs() < 1e-6);
        assert_eq!(s.reinforcing, 0.0);
        assert_eq!(s.feedback, 0.0);
        assert_eq!(s.lateral, 0.0);
        assert!((s.bystander - 0.1).abs() < 1e-6);
    }

    /// `transition_window = 1` extends transition_like to adjacent bins.
    #[test]
    fn compute_role_scores_one_transition_window_adjacent_bins() {
        let bif = BifurcationResult {
            b_star: Some(4), // 0-based bin 4 ↔ 1-based bin 5
            reason: BifurcationReason::Crossing,
        };
        let params = RoleScoreParams {
            tau_0: 1.0,
            h_fork: 0.5,
            min_n_eff_cells: 0.0,
            transition_window: 1,
        };
        let at_star = compute_role_scores_one(&bin_stub(5, 0), 0.8, 0.2, 1.0, bif, &params);
        assert!((at_star.transition_like - 0.8).abs() < 1e-6);
        let adjacent = compute_role_scores_one(&bin_stub(4, 0), 0.8, 0.2, 1.0, bif, &params);
        assert!((adjacent.transition_like - 0.8).abs() < 1e-6);
        let outside = compute_role_scores_one(&bin_stub(2, 0), 0.8, 0.2, 1.0, bif, &params);
        assert_eq!(outside.transition_like, 0.0);
    }

    /// `lateral` is the only score that consumes `spec_max_other`; for
    /// a pre-fork row with Δτ << 0 (so σ(-Δτ/τ_0) ≈ 1), lateral should
    /// approach `spec_max_other` while instructive stays at 0 because
    /// `spec_self = 0`. This tests the weak-tier routing.
    #[test]
    fn compute_role_scores_one_lateral_uses_opposing_spec() {
        let bin = bin_stub(/*b=*/ 1, 0);
        let bif = BifurcationResult {
            b_star: Some(4),
            reason: BifurcationReason::Crossing,
        };
        let params = RoleScoreParams {
            tau_0: 1.0,
            h_fork: 0.5,
            min_n_eff_cells: 0.0,
            transition_window: 0,
        };
        let s = compute_role_scores_one(
            &bin, /*spec_self*/ 0.0, // isolate lateral
            /*spec_max_other*/ 0.7,
            /*delta_tau*/ -5.0, // σ(-(-5)) = σ(5) ≈ 0.9933
            bif, &params,
        );
        assert_eq!(s.instructive, 0.0);
        assert_eq!(s.reinforcing, 0.0);
        assert_eq!(s.feedback, 0.0);
        // lateral = 0.7 · 1 · σ(5) ≈ 0.7 · 0.9933 ≈ 0.6953
        assert!(
            (s.lateral - 0.7 * 0.9933).abs() < 1e-3,
            "lateral = {}",
            s.lateral
        );
        assert!((s.bystander - (1.0 - s.lateral)).abs() < 1e-6);
    }

    /// `bystander = 1 - max(instructive, transition_like, reinforcing, feedback, lateral)`
    /// algebraically; verify on a row that fires multiple main terms at
    /// different magnitudes so the `max` picks correctly.
    #[test]
    fn compute_role_scores_one_bystander_is_one_minus_max_main() {
        let bin = bin_stub(/*b=*/ 5, 0);
        let bif = BifurcationResult {
            b_star: Some(2),
            reason: BifurcationReason::Crossing,
        };
        let params = RoleScoreParams {
            tau_0: 1.0,
            h_fork: 0.5,
            min_n_eff_cells: 0.0,
            transition_window: 0,
        };
        // b > b* so only post-fork terms fire.
        // reinforcing = 0.9 · 1 · σ(2) ≈ 0.9 · 0.8808 ≈ 0.7927
        // feedback    = 0.9 · 1 · σ(-2) ≈ 0.9 · 0.1192 ≈ 0.1073
        let s = compute_role_scores_one(
            &bin, /*spec_self*/ 0.9, /*spec_max_other*/ 0.9, /*delta_tau*/ 2.0, bif,
            &params,
        );
        let expected_max = s
            .reinforcing
            .max(s.instructive)
            .max(s.feedback)
            .max(s.lateral);
        assert!(
            (s.bystander - (1.0 - expected_max)).abs() < 1e-6,
            "bystander {} != 1 - max {}",
            s.bystander,
            expected_max
        );
        assert!(
            (expected_max - s.reinforcing).abs() < 1e-6,
            "max should be reinforcing, got {}",
            expected_max
        );
    }

    /// Sigmoid smoke test on the algorithmic boundary: σ(0) = 0.5
    /// exactly, σ(large positive) ≈ 1, σ(large negative) ≈ 0, σ(NaN)
    /// propagates NaN.  Indirect — exercised via the zero-Δτ and
    /// very-positive-Δτ scoring calls — but worth an explicit row.
    #[test]
    fn compute_role_scores_one_nan_delta_tau_propagates() {
        let bin = bin_stub(/*b=*/ 5, 0);
        let bif = BifurcationResult {
            b_star: Some(2),
            reason: BifurcationReason::Crossing,
        };
        let params = RoleScoreParams {
            tau_0: 1.0,
            h_fork: 0.5,
            min_n_eff_cells: 0.0,
            transition_window: 0,
        };
        let s = compute_role_scores_one(
            &bin,
            /*spec_self*/ 0.5,
            /*spec_max_other*/ 0.5,
            /*delta_tau*/ f32::NAN,
            bif,
            &params,
        );
        // Only post-fork terms would normally fire (b > b*).  NaN Δτ
        // flows through σ(NaN / 1) = NaN, so both reinforcing and
        // feedback should be NaN on this row.
        assert!(s.reinforcing.is_nan(), "reinforcing = {}", s.reinforcing);
        assert!(s.feedback.is_nan(), "feedback    = {}", s.feedback);
        // instructive / lateral would multiply by pre_ind = 0 — however
        // 0 * NaN = NaN in IEEE-754, so they also propagate NaN.  This
        // documents the observed behavior: callers must filter NaN Δτ
        // rows before running the taxonomy.
        assert!(s.instructive.is_nan());
        assert!(s.lateral.is_nan());
        // bystander = 1 - max(all NaN) — max of NaNs is NaN, so
        // bystander is also NaN.
        assert!(s.bystander.is_nan());
    }

    // ─── Row-level n_eff_cells gate tests ───────────────────────────────────

    /// Helper extending `bin_stub` with an explicit `n_eff_cells` value
    /// Role scoring tests set `min_n_eff_cells = 0.0`; gate-specific tests vary
    /// this field explicitly.
    fn bin_stub_with_n_eff(bin_idx: usize, lineage_idx: usize, n_eff_cells: f32) -> LrBinResult {
        let mut b = bin_stub(bin_idx, lineage_idx);
        b.n_eff_cells = n_eff_cells;
        b
    }

    /// A below-threshold row (`n_eff_cells = 4.0 < min_n_eff_cells = 5.0`)
    /// collapses to `(0,0,0,0,0,1)` even with a high-scoring configuration.
    /// The reliability gate is applied before formula evaluation.
    #[test]
    fn compute_role_scores_one_p3_gate_zeroes_unreliable_row() {
        let bin = bin_stub_with_n_eff(/*b=*/ 2, 0, /*n_eff*/ 4.0);
        let bif = BifurcationResult {
            b_star: Some(5), // b < b*
            reason: BifurcationReason::Crossing,
        };
        let params = RoleScoreParams {
            tau_0: 1.0,
            h_fork: 0.5,
            min_n_eff_cells: 5.0,
            transition_window: 0,
        };
        let s = compute_role_scores_one(
            &bin, /*spec_self*/ 0.9, /*spec_max_other*/ 0.9, /*delta_tau*/ 5.0, bif,
            &params,
        );
        assert_eq!(s.instructive, 0.0, "gate must zero instructive");
        assert_eq!(s.reinforcing, 0.0);
        assert_eq!(s.feedback, 0.0);
        assert_eq!(s.lateral, 0.0);
        assert_eq!(s.bystander, 1.0, "gate must set bystander = 1");
    }

    /// At-or-above-threshold row: `n_eff_cells = 6.0 >= 5.0` passes the
    /// gate, so the same hot scenario as the unreliable test above
    /// produces normal  scoring.  Pre-fork forward → instructive
    /// ≈ 0.9 · σ(5) ≈ 0.8940.  This is the matched-pair counterpart to
    /// `_p3_gate_zeroes_unreliable_row` — same inputs except n_eff.
    #[test]
    fn compute_role_scores_one_p3_gate_passes_reliable_row() {
        let bin = bin_stub_with_n_eff(/*b=*/ 2, 0, /*n_eff*/ 6.0);
        let bif = BifurcationResult {
            b_star: Some(5),
            reason: BifurcationReason::Crossing,
        };
        let params = RoleScoreParams {
            tau_0: 1.0,
            h_fork: 0.5,
            min_n_eff_cells: 5.0,
            transition_window: 0,
        };
        let s = compute_role_scores_one(
            &bin, /*spec_self*/ 0.9, /*spec_max_other*/ 0.0, // isolate instructive
            /*delta_tau*/ 5.0, bif, &params,
        );
        // instructive = 0.9 · 1 · σ(5) ≈ 0.9 · 0.9933 ≈ 0.8940
        assert!(
            (s.instructive - 0.9 * 0.9933).abs() < 1e-3,
            "instructive = {} (expected ≈ 0.8940)",
            s.instructive
        );
        assert_eq!(s.reinforcing, 0.0);
        assert_eq!(s.feedback, 0.0);
        assert_eq!(s.lateral, 0.0); // spec_max_other = 0
        assert!(
            s.bystander < 0.2,
            "bystander = {} should be small",
            s.bystander
        );
    }

    /// Boundary semantics: `min_n_eff_cells = 0.0` disables the gate
    /// entirely.  A row with `n_eff_cells = 0.0` (the worst possible)
    /// MUST then score normally under the role-scoring equations. This locks the documented
    /// disable-path for unit-test isolation and for small datasets where
    /// a 5.0 minimum would be too strict.
    #[test]
    fn compute_role_scores_one_p3_gate_disabled_at_zero_threshold() {
        let bin = bin_stub_with_n_eff(/*b=*/ 2, 0, /*n_eff*/ 0.0);
        let bif = BifurcationResult {
            b_star: Some(5),
            reason: BifurcationReason::Crossing,
        };
        let params = RoleScoreParams {
            tau_0: 1.0,
            h_fork: 0.5,
            min_n_eff_cells: 0.0,
            transition_window: 0, // disabled
        };
        let s = compute_role_scores_one(
            &bin, /*spec_self*/ 0.9, /*spec_max_other*/ 0.0, /*delta_tau*/ 5.0, bif,
            &params,
        );
        // Same expectation as the reliable-row test: gate off ⇒
        //  fires normally regardless of n_eff_cells value.
        assert!(
            (s.instructive - 0.9 * 0.9933).abs() < 1e-3,
            "disabled but instructive = {}",
            s.instructive
        );
    }
}
