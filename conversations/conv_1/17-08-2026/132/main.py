#!/usr/bin/env python3

import math
from collections import defaultdict

# =============================================================================
# EXPERIMENT 197
# EXACT TERM-LEVEL k-DEPENDENCE EXTRACTION
# =============================================================================
#
# Goal:
#   Determine whether the k-dependence is already absent from the current
#   exact pre-convolution summands, or whether it can enter through:
#
#       1. j-range / truncation
#       2. binomial indices
#       3. multiplicative factors
#       4. signs
#
# This is deliberately one self-contained main.py.
#
# IMPORTANT:
#   We do NOT assume that the ordinary convolution is the final pq formula.
#   We use the known exact anchors only as external diagnostics.
#
# =============================================================================


# -----------------------------------------------------------------------------
# BASIC UTILITIES
# -----------------------------------------------------------------------------

def C(n, r):
    if n < 0 or r < 0 or r > n:
        return 0
    return math.comb(n, r)


def sign_pm1(x):
    return -1 if x % 2 else 1


# -----------------------------------------------------------------------------
# KNOWN EXACT ANCHORS FROM EXPERIMENTS 191-194
# -----------------------------------------------------------------------------

ANCHORS = [
    ("A1", 1, 10, 2, 27),
    ("A2", 3, 20, 6, 935),
    ("A3", 7, 40, 15, -1797818),
    ("B1", 3, 7, 3, -3),
    ("B2", 3, 9, 4, 9),
    ("B3", 3, 11, 5, -16),
    ("B4", 5, 11, 5, -5),
    ("B5", 5, 19, 9, -196),
    ("B6", 9, 21, 11, -84),
]


# -----------------------------------------------------------------------------
# REFERENCE CONVOLUTION
# -----------------------------------------------------------------------------

def reference_terms(k, ell, s):
    """
    The convolution already identified in Experiments 185-186:

        sum_j (-1)^j C(ell,j) C(ell-j,s-j)

    k is deliberately unused.

    This is exactly the point of the experiment:
    measure where the loss of k occurs.
    """

    out = []

    for j in range(0, ell + 1):

        a = C(ell, j)
        b = C(ell - j, s - j)

        if a == 0 or b == 0:
            continue

        sg = sign_pm1(j)
        value = sg * a * b

        out.append({
            "j": j,
            "sign": sg,
            "A": a,
            "B": b,
            "value": value,

            "A_top": ell,
            "A_bot": j,

            "B_top": ell - j,
            "B_bot": s - j,
        })

    return out


def term_sum(terms):
    return sum(t["value"] for t in terms)


# -----------------------------------------------------------------------------
# EXPERIMENTAL PRE-COLLAPSE STRUCTURES
# -----------------------------------------------------------------------------
#
# We now explicitly construct the finite structural variants already isolated
# in Experiments 187-188.
#
# These are NOT claimed to be the pq theorem.
#
# They are diagnostic mechanisms for locating where k can enter.
#
# -----------------------------------------------------------------------------

def variant_A(k, ell, s):
    """
    Second upper parameter + k.
    """
    out = []

    for j in range(0, ell + 1):
        a = C(ell, j)
        b = C(ell - j + k, s - j)

        if a == 0 or b == 0:
            continue

        sg = sign_pm1(j)

        out.append({
            "j": j,
            "sign": sg,
            "A": a,
            "B": b,
            "value": sg * a * b,
            "A_top": ell,
            "A_bot": j,
            "B_top": ell - j + k,
            "B_bot": s - j,
        })

    return out


def variant_B(k, ell, s):
    """
    Target index shifted by -k.
    """
    out = []

    target = s - k

    for j in range(0, ell + 1):
        a = C(ell, j)
        b = C(ell - j, target - j)

        if a == 0 or b == 0:
            continue

        sg = sign_pm1(j)

        out.append({
            "j": j,
            "sign": sg,
            "A": a,
            "B": b,
            "value": sg * a * b,
            "A_top": ell,
            "A_bot": j,
            "B_top": ell - j,
            "B_bot": target - j,
        })

    return out


def variant_C(k, ell, s):
    """
    Truncate the j-range at j <= k.
    """
    out = []

    for j in range(0, min(ell, k) + 1):
        a = C(ell, j)
        b = C(ell - j, s - j)

        if a == 0 or b == 0:
            continue

        sg = sign_pm1(j)

        out.append({
            "j": j,
            "sign": sg,
            "A": a,
            "B": b,
            "value": sg * a * b,
            "A_top": ell,
            "A_bot": j,
            "B_top": ell - j,
            "B_bot": s - j,
        })

    return out


def variant_D(k, ell, s):
    """
    Second upper +k AND target s-k.
    """
    out = []

    target = s - k

    for j in range(0, ell + 1):
        a = C(ell, j)
        b = C(ell - j + k, target - j)

        if a == 0 or b == 0:
            continue

        sg = sign_pm1(j)

        out.append({
            "j": j,
            "sign": sg,
            "A": a,
            "B": b,
            "value": sg * a * b,
            "A_top": ell,
            "A_bot": j,
            "B_top": ell - j + k,
            "B_bot": target - j,
        })

    return out


def variant_E(k, ell, s):
    """
    First upper -k.
    """
    out = []

    for j in range(0, ell + 1):
        a = C(ell - k, j)
        b = C(ell - j, s - j)

        if a == 0 or b == 0:
            continue

        sg = sign_pm1(j)

        out.append({
            "j": j,
            "sign": sg,
            "A": a,
            "B": b,
            "value": sg * a * b,
            "A_top": ell - k,
            "A_bot": j,
            "B_top": ell - j,
            "B_bot": s - j,
        })

    return out


def variant_F(k, ell, s):
    """
    First upper +k.
    """
    out = []

    for j in range(0, ell + 1):
        a = C(ell + k, j)
        b = C(ell - j, s - j)

        if a == 0 or b == 0:
            continue

        sg = sign_pm1(j)

        out.append({
            "j": j,
            "sign": sg,
            "A": a,
            "B": b,
            "value": sg * a * b,
            "A_top": ell + k,
            "A_bot": j,
            "B_top": ell - j,
            "B_bot": s - j,
        })

    return out


VARIANTS = {
    "reference": reference_terms,
    "A_second_upper_plus_k": variant_A,
    "B_target_minus_k": variant_B,
    "C_truncate_j_at_k": variant_C,
    "D_both": variant_D,
    "E_first_upper_minus_k": variant_E,
    "F_first_upper_plus_k": variant_F,
}


# -----------------------------------------------------------------------------
# ANCHOR MATCHING
# -----------------------------------------------------------------------------

def evaluate_anchor_variant(fn, k, ell, s):
    return term_sum(fn(k, ell, s))


def anchor_matches(fn):
    matches = 0

    for _, k, ell, s, expected in ANCHORS:
        value = evaluate_anchor_variant(fn, k, ell, s)
        if value == expected:
            matches += 1

    return matches


# -----------------------------------------------------------------------------
# TERM-LEVEL k COMPARISON
# -----------------------------------------------------------------------------

def compare_fixed_group(fn, ell, s, ks):

    data = {}

    for k in ks:
        terms = fn(k, ell, s)

        by_j = {}
        for t in terms:
            by_j[t["j"]] = t

        data[k] = by_j

    all_j = sorted(set(
        j
        for by_j in data.values()
        for j in by_j
    ))

    diagnostics = {
        "support": False,
        "sign": False,
        "A_factor": False,
        "B_factor": False,
        "index": False,
        "value": False,
    }

    for j in all_j:

        entries = []

        for k in ks:
            t = data[k].get(j)

            if t is None:
                entries.append(None)
            else:
                entries.append(t)

        # Support / truncation.
        if any(t is None for t in entries):
            diagnostics["support"] = True

        present = [t for t in entries if t is not None]

        if not present:
            continue

        # Sign.
        signs = set(t["sign"] for t in present)
        if len(signs) > 1:
            diagnostics["sign"] = True

        # Factors.
        As = set(t["A"] for t in present)
        Bs = set(t["B"] for t in present)

        if len(As) > 1:
            diagnostics["A_factor"] = True

        if len(Bs) > 1:
            diagnostics["B_factor"] = True

        # Binomial index structure.
        idx = set(
            (
                t["A_top"],
                t["A_bot"],
                t["B_top"],
                t["B_bot"],
            )
            for t in present
        )

        if len(idx) > 1:
            diagnostics["index"] = True

        # Actual term values.
        vals = set(t["value"] for t in present)
        if len(vals) > 1:
            diagnostics["value"] = True

    return diagnostics, data


# -----------------------------------------------------------------------------
# INTERESTING k-NONZERO STRUCTURE SEARCH
# -----------------------------------------------------------------------------
#
# We search ONLY small local shifts.
# This is intentionally bounded, unlike Experiments 192-194.
#
# Structure:
#
#   C(ell + dk, j + dj) C(ell - j + dk2, s - j + ds)
#
# with very small shifts.
# -----------------------------------------------------------------------------

def shifted_family(k, ell, s, dk1, dj1, dk2, ds2):
    out = []

    for j in range(0, ell + 1):

        a = C(ell + dk1 * k, j + dj1 * k)

        b = C(
            ell - j + dk2 * k,
            s - j + ds2 * k
        )

        if a == 0 or b == 0:
            continue

        sg = sign_pm1(j)

        out.append({
            "j": j,
            "sign": sg,
            "A": a,
            "B": b,
            "value": sg * a * b,

            "A_top": ell + dk1 * k,
            "A_bot": j + dj1 * k,

            "B_top": ell - j + dk2 * k,
            "B_bot": s - j + ds2 * k,
        })

    return out


# -----------------------------------------------------------------------------
# SMALL SHIFT SEARCH
# -----------------------------------------------------------------------------

def small_shift_search():

    print()
    print("=" * 78)
    print("6. BOUNDED LOCAL PRE-CONVOLUTION SHIFT SEARCH")
    print("=" * 78)

    candidates = []

    # Very small range only.
    for dk1 in [-1, 0, 1]:
        for dj1 in [-1, 0, 1]:
            for dk2 in [-1, 0, 1]:
                for ds2 in [-1, 0, 1]:

                    def fn(k, ell, s,
                           a=dk1, b=dj1, c=dk2, d=ds2):
                        return shifted_family(k, ell, s, a, b, c, d)

                    matches = anchor_matches(fn)

                    # Keep only candidates that explain at least one anchor.
                    if matches > 0:

                        candidates.append(
                            (
                                matches,
                                dk1,
                                dj1,
                                dk2,
                                ds2,
                            )
                        )

    candidates.sort(reverse=True)

    for matches, dk1, dj1, dk2, ds2 in candidates[:20]:

        print(
            f"dk1={dk1:+d} dj1={dj1:+d} "
            f"dk2={dk2:+d} ds2={ds2:+d}"
            f"    anchors={matches}/9"
        )

    if not candidates:
        print("No local shift reproduced even one anchor.")

    return candidates


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 197")
    print("EXACT TERM-LEVEL k-DEPENDENCE EXTRACTION")
    print("=" * 78)

    # -------------------------------------------------------------------------
    # 1. REFERENCE CONVOLUTION
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. REFERENCE CONVOLUTION")
    print("=" * 78)

    tests = [
        (1, 5, 2),
        (1, 6, 3),
        (3, 8, 3),
        (3, 9, 4),
        (5, 10, 5),
    ]

    for k, ell, s in tests:

        terms = reference_terms(k, ell, s)

        print(
            f"k={k:2d} ell={ell:2d} s={s:2d}"
            f" sum={term_sum(terms):8d}"
        )

    # -------------------------------------------------------------------------
    # 2. ANCHOR AUDIT
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. CURRENT STRUCTURAL VARIANTS AGAINST 9 ANCHORS")
    print("=" * 78)

    for name, fn in VARIANTS.items():

        matches = anchor_matches(fn)

        print(
            f"{name:30s}"
            f" {matches}/9"
        )

    # -------------------------------------------------------------------------
    # 3. FIXED ell,s VARY k
    # -------------------------------------------------------------------------

    groups = [
        (8, 3, [1, 3, 5, 7]),
        (8, 4, [1, 3, 5, 7]),
        (9, 4, [1, 3, 5, 7, 9]),
        (10, 4, [1, 3, 5, 7, 9]),
        (10, 5, [1, 3, 5, 7, 9]),
        (12, 5, [1, 3, 5, 7, 9, 11]),
        (12, 6, [1, 3, 5, 7, 9, 11]),
        (20, 6, [1, 3, 5, 7, 9, 11]),
        (20, 7, [1, 3, 5, 7, 9, 11]),
        (40, 15, [1, 3, 5, 7, 9, 11]),
    ]

    print()
    print("=" * 78)
    print("3. WHERE DOES k ENTER IN EACH STRUCTURE?")
    print("=" * 78)

    for name, fn in VARIANTS.items():

        print()
        print("-" * 78)
        print(name)

        totals = defaultdict(int)

        for ell, s, ks in groups:

            diag, _ = compare_fixed_group(
                fn,
                ell,
                s,
                ks,
            )

            for key, value in diag.items():
                if value:
                    totals[key] += 1

        print(
            f"  support-varying groups = "
            f"{totals['support']}/{len(groups)}"
        )

        print(
            f"  sign-varying groups    = "
            f"{totals['sign']}/{len(groups)}"
        )

        print(
            f"  A-factor-varying       = "
            f"{totals['A_factor']}/{len(groups)}"
        )

        print(
            f"  B-factor-varying       = "
            f"{totals['B_factor']}/{len(groups)}"
        )

        print(
            f"  index-varying          = "
            f"{totals['index']}/{len(groups)}"
        )

        print(
            f"  term-value-varying     = "
            f"{totals['value']}/{len(groups)}"
        )

    # -------------------------------------------------------------------------
    # 4. EXPLICIT TERM COMPARISON
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXPLICIT TERM COMPARISON FOR DIAGNOSTIC CASE")
    print("=" * 78)

    ell = 12
    s = 5
    ks = [1, 3, 5, 7, 9]

    for name, fn in VARIANTS.items():

        print()
        print(name)

        for k in ks:

            terms = fn(k, ell, s)

            print(
                f"  k={k:2d}"
                f" total={term_sum(terms):10d}"
                f" terms={len(terms):2d}"
            )

        diag, data = compare_fixed_group(
            fn,
            ell,
            s,
            ks,
        )

        print(f"  support varies = {diag['support']}")
        print(f"  sign varies    = {diag['sign']}")
        print(f"  A factors vary = {diag['A_factor']}")
        print(f"  B factors vary = {diag['B_factor']}")
        print(f"  indices vary   = {diag['index']}")
        print(f"  values vary    = {diag['value']}")

    # -------------------------------------------------------------------------
    # 5. ANCHOR TERM PATTERN
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. ANCHOR VALUE / REFERENCE VALUE DIAGNOSTIC")
    print("=" * 78)

    for name, k, ell, s, exact in ANCHORS:

        ref = term_sum(
            reference_terms(k, ell, s)
        )

        print(
            f"{name:3s}"
            f" k={k:2d}"
            f" ell={ell:2d}"
            f" s={s:2d}"
            f" exact={exact:12d}"
            f" ref={ref:8d}"
        )

    # -------------------------------------------------------------------------
    # 6. BOUNDED LOCAL SHIFT SEARCH
    # -------------------------------------------------------------------------

    candidates = small_shift_search()

    # -------------------------------------------------------------------------
    # 7. FINAL DIAGNOSTIC
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    if candidates:

        best = candidates[0]

        print(
            "A local pre-convolution structure produced nonzero anchor agreement."
        )

        print(
            "Best local candidate:"
            f" anchors={best[0]}/9,"
            f" dk1={best[1]:+d},"
            f" dj1={best[2]:+d},"
            f" dk2={best[3]:+d},"
            f" ds2={best[4]:+d}"
        )

        print()
        print(
            "IMPORTANT: this is only a diagnostic survivor."
        )

        print(
            "It must be derived from the exact pq construction before"
            " it can be accepted as mathematical evidence."
        )

    else:

        print(
            "No bounded pre-convolution shift reproduced an anchor."
        )

        print()
        print(
            "This strengthens the conclusion from Experiments 187-196:"
        )

        print(
            "the missing k-dependence is probably not a small affine"
            " perturbation of the current convolution."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "derive the exact pq-level combinatorial object before"
            " the binomial convolution is introduced."
        )

    print()
    print(
        "The critical distinction remains:"
    )

    print(
        "    pq structure -> exact summand -> admissible j-range"
    )

    print(
        "not"
    )

    print(
        "    reduced (k,ell,s) state -> guessed convolution."
    )

    print("=" * 78)


if __name__ == "__main__":
    main()