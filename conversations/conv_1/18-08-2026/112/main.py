#!/usr/bin/env python3

from fractions import Fraction
from math import factorial, gcd
from functools import reduce


# ============================================================================
# EXPERIMENT 115 — EXACT CORRECTED SECOND FALLING-FACTORIAL AUDIT
# ============================================================================
#
# Experiment 114 found:
#
#     C_p(k) = (K-k)_s Q_p(k)
#
# but its second falling-factorial reconstruction failed because the
# transform accidentally used Python floor division on Fraction objects.
#
# This experiment:
#
#   1. reconstructs every Q_p(k) correctly in the k_(r) basis;
#   2. independently verifies the Newton/falling transform;
#   3. checks exact triangular support;
#   4. inspects the transformed coefficient matrices;
#   5. searches for terminal/lower-end factorial factors;
#   6. compares A/B sectors only after the transform is proven correct.
#
# ALL arithmetic is exact QQ.
# No SymPy.
# No floating point.
#
# ============================================================================


# ----------------------------------------------------------------------------
# DATA
# ----------------------------------------------------------------------------

A_even = {
    0: [Fraction(-12879), Fraction(-28241), Fraction(-26989),
        Fraction(-13611), Fraction(-17875, 6), Fraction(-116923, 1680),
        Fraction(-5, 144)],
    2: [Fraction(2797337, 1920), Fraction(8986567, 2880),
        Fraction(36298273, 13440), Fraction(11670379, 11520),
        Fraction(19954213, 161280), Fraction(-9389, 4032), Fraction(0)],
    4: [Fraction(-2083937, 30720), Fraction(-83529, 640),
        Fraction(-7965025, 129024), Fraction(2225141, 46080),
        Fraction(710501, 215040), Fraction(0), Fraction(0)],
    6: [Fraction(85591, 61440), Fraction(83651, 46080),
        Fraction(-3174439, 2580480), Fraction(0), Fraction(0),
        Fraction(0), Fraction(0)],
    8: [Fraction(-4913, 491520), Fraction(0), Fraction(0),
        Fraction(0), Fraction(0), Fraction(0), Fraction(0)],
}

A_odd = {
    1: [Fraction(12143, 560), Fraction(2699231, 1344),
        Fraction(9120441, 2240), Fraction(21975383, 6720),
        Fraction(4460869, 5760), Fraction(42929, 1680), Fraction(0)],
    3: [Fraction(-989, 11520), Fraction(-12024227, 46080),
        Fraction(-175956721, 322560), Fraction(-67903883, 161280),
        Fraction(-2590159, 53760), Fraction(0), Fraction(0)],
    5: [Fraction(-517, 23040), Fraction(396119, 30720),
        Fraction(26625517, 1290240), Fraction(-234707, 129024),
        Fraction(0), Fraction(0), Fraction(0)],
    7: [Fraction(373, 1290240), Fraction(-1028053, 5160960),
        Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(0)],
}

B_even = {
    0: [Fraction(12980463, 1024), Fraction(6255583, 256),
        Fraction(87841139, 4480), Fraction(16998339, 2240),
        Fraction(1091983, 896), Fraction(1553, 240)],
    2: [Fraction(-19344659, 15360), Fraction(-26986999, 11520),
        Fraction(-2066529, 1120), Fraction(-573325, 576),
        Fraction(3312053, 40320), Fraction(0)],
    4: [Fraction(129415, 3072), Fraction(267779, 3840),
        Fraction(224417, 4480), Fraction(-101119, 5040),
        Fraction(0), Fraction(0)],
    6: [Fraction(-2267, 5120), Fraction(-6053, 11520),
        Fraction(0), Fraction(0), Fraction(0), Fraction(0)],
}

B_odd = {
    1: [Fraction(-584531, 35840), Fraction(-2908483, 1920),
        Fraction(-31233169, 13440), Fraction(-1446167, 1344),
        Fraction(-22259149, 40320), Fraction(-301, 240)],
    3: [Fraction(-59257, 46080), Fraction(186547, 1440),
        Fraction(367433, 1680), Fraction(126549, 448),
        Fraction(-162139, 40320), Fraction(0)],
    5: [Fraction(4457, 46080), Fraction(-16819, 5760),
        Fraction(-5769, 896), Fraction(0), Fraction(0), Fraction(0)],
    7: [Fraction(-421, 322560), Fraction(0), Fraction(0),
        Fraction(0), Fraction(0), Fraction(0)],
}


CHANNELS = {
    "A-even": {"data": A_even, "K": 6, "c": 0},
    "A-odd": {"data": A_odd, "K": 6, "c": 0},
    "B-even": {"data": B_even, "K": 5, "c": 1},
    "B-odd": {"data": B_odd, "K": 5, "c": 1},
}


# ----------------------------------------------------------------------------
# POLYNOMIAL UTILITIES
# ----------------------------------------------------------------------------

def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p if p else [Fraction(0)]


def zero_poly():
    return [Fraction(0)]


def is_zero(p):
    return all(c == 0 for c in p)


def degree(p):
    p = trim(p)
    return -1 if is_zero(p) else len(p) - 1


def add(a, b):
    n = max(len(a), len(b))
    out = [Fraction(0)] * n
    for i in range(n):
        if i < len(a):
            out[i] += a[i]
        if i < len(b):
            out[i] += b[i]
    return trim(out)


def sub(a, b):
    n = max(len(a), len(b))
    out = [Fraction(0)] * n
    for i in range(n):
        if i < len(a):
            out[i] += a[i]
        if i < len(b):
            out[i] -= b[i]
    return trim(out)


def scale(a, q):
    return trim([Fraction(q) * x for x in a])


def mul(a, b):
    if is_zero(a) or is_zero(b):
        return zero_poly()

    out = [Fraction(0)] * (len(a) + len(b) - 1)

    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj

    return trim(out)


def eval_poly(p, x):
    x = Fraction(x)
    acc = Fraction(0)

    for c in reversed(p):
        acc = acc * x + c

    return acc


def divmod_exact(f, g):
    f = trim(f)
    g = trim(g)

    if is_zero(g):
        raise ZeroDivisionError("zero polynomial divisor")

    df = degree(f)
    dg = degree(g)

    if df < dg:
        return zero_poly(), f[:]

    q = [Fraction(0)] * (df - dg + 1)
    r = f[:]

    while not is_zero(r) and degree(r) >= dg:
        dr = degree(r)
        coeff = r[-1] / g[-1]
        shift = dr - dg

        q[shift] += coeff

        term = [Fraction(0)] * shift + [coeff]
        r = sub(r, mul(term, g))

    return trim(q), trim(r)


def div_exact(f, g):
    q, r = divmod_exact(f, g)

    if not is_zero(r):
        raise ValueError("non-exact polynomial division")

    return q


def gcd_poly(a, b):
    a = trim(a)
    b = trim(b)

    while not is_zero(b):
        _, r = divmod_exact(a, b)
        a, b = b, r

    if is_zero(a):
        return [Fraction(0)]

    lead = a[-1]
    return scale(a, Fraction(1, 1) / lead)


# ----------------------------------------------------------------------------
# INTERPOLATION
# ----------------------------------------------------------------------------

def interpolate(values):
    K = len(values) - 1
    result = zero_poly()

    for i, yi in enumerate(values):
        basis = [Fraction(1)]
        denom = Fraction(1)

        for j in range(K + 1):
            if j == i:
                continue

            basis = mul(basis, [Fraction(-j), Fraction(1)])
            denom *= Fraction(i - j)

        result = add(result, scale(basis, yi / denom))

    return trim(result)


# ----------------------------------------------------------------------------
# FALLING FACTORIAL POLYNOMIALS
# ----------------------------------------------------------------------------

def falling(r):
    p = [Fraction(1)]

    for j in range(r):
        p = mul(p, [Fraction(-j), Fraction(1)])

    return p


def terminal(K, s):
    p = [Fraction(1)]

    for j in range(s):
        # K-k-j
        p = mul(p, [Fraction(K - j), Fraction(-1)])

    return p


# ----------------------------------------------------------------------------
# SUPPORT LAW
# ----------------------------------------------------------------------------

def support_s(p, c):
    return max(
        p - 2,
        (p - c + 1) // 2,
    )


# ----------------------------------------------------------------------------
# CORRECT FALLING-FACTORIAL TRANSFORM
#
# If
#
#     Q(k) = q_0 + q_1 k_(1) + ... + q_d k_(d),
#
# then at k=n:
#
#     Q(n) = sum_{r=0}^n q_r * n!/(n-r)!.
#
# IMPORTANT:
#
#     factorial(n) / factorial(n-r)
#
# is represented as Fraction(...), NOT floor division.
# ----------------------------------------------------------------------------

def to_falling_basis(Q):
    d = degree(Q)

    if d < 0:
        return []

    coeffs = [Fraction(0)] * (d + 1)

    for n in range(d + 1):
        residual = eval_poly(Q, n)

        for r in range(n):
            falling_value = Fraction(
                factorial(n),
                factorial(n - r),
            )

            residual -= coeffs[r] * falling_value

        coeffs[n] = residual / factorial(n)

    return trim(coeffs)


def from_falling_basis(coeffs):
    result = zero_poly()

    for r, q in enumerate(coeffs):
        if q:
            result = add(
                result,
                scale(falling(r), q),
            )

    return trim(result)


# ----------------------------------------------------------------------------
# INTEGER SIGNATURE
# ----------------------------------------------------------------------------

def primitive_signature(values):
    values = trim(values)

    if is_zero(values):
        return [0]

    den_lcm = 1

    for x in values:
        den = x.denominator
        den_lcm = den_lcm * den // gcd(den_lcm, den)

    ints = [int(x * den_lcm) for x in values]

    g = reduce(
        gcd,
        [abs(v) for v in ints if v != 0],
        0,
    )

    if g:
        ints = [v // g for v in ints]

    # Canonical sign: first nonzero positive.
    for v in ints:
        if v:
            if v < 0:
                ints = [-x for x in ints]
            break

    return ints


# ----------------------------------------------------------------------------
# MATRIX RANK
# ----------------------------------------------------------------------------

def rank(matrix):
    if not matrix:
        return 0

    A = [
        [Fraction(x) for x in row]
        for row in matrix
    ]

    width = max(len(row) for row in A)

    for row in A:
        row.extend([Fraction(0)] * (width - len(row)))

    m = len(A)
    n = width

    r = 0

    for c in range(n):
        pivot = None

        for i in range(r, m):
            if A[i][c] != 0:
                pivot = i
                break

        if pivot is None:
            continue

        A[r], A[pivot] = A[pivot], A[r]

        p = A[r][c]
        A[r] = [x / p for x in A[r]]

        for i in range(m):
            if i == r:
                continue

            f = A[i][c]

            if f:
                A[i] = [
                    A[i][j] - f * A[r][j]
                    for j in range(n)
                ]

        r += 1

        if r == m:
            break

    return r


# ----------------------------------------------------------------------------
# GCD / FACTOR SEARCH
# ----------------------------------------------------------------------------

def lower_factor(r):
    return falling(r)


def shifted_factor(a, r):
    p = [Fraction(1)]

    for j in range(r):
        p = mul(
            p,
            [Fraction(-(a + j)), Fraction(1)]
        )

    return p


def exact_factor_hits(Q):
    d = degree(Q)

    if d < 1:
        return []

    hits = []

    for r in range(1, d + 1):
        f = lower_factor(r)
        _, rem = divmod_exact(Q, f)

        if is_zero(rem):
            hits.append(("k_(r)", r))

    for a in range(-3, 4):
        for r in range(1, d + 1):
            f = shifted_factor(a, r)
            _, rem = divmod_exact(Q, f)

            if is_zero(rem):
                hits.append((f"(k-{a})_(r)", r))

    return hits


# ----------------------------------------------------------------------------
# PRECOMPUTED TEST TRANSFORM
#
# This independently verifies the transform on the known basis itself:
#
#     k_(r) -> coefficient vector with a single 1 at r.
# ----------------------------------------------------------------------------

def transform_self_test():
    for d in range(0, 8):
        for r in range(d + 1):
            basis_poly = falling(r)
            coeffs = to_falling_basis(basis_poly)

            expected = [Fraction(0)] * (d + 1)
            expected[r] = Fraction(1)
            expected = trim(expected)

            if coeffs != expected:
                return False, d, r, coeffs, expected

            reconstructed = from_falling_basis(coeffs)

            if reconstructed != trim(basis_poly):
                return False, d, r, reconstructed, basis_poly

    return True, None, None, None, None


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 78)
    print("EXPERIMENT 115 — EXACT CORRECTED SECOND FALLING-FACTORIAL AUDIT")
    print("=" * 78)
    print()

    # ------------------------------------------------------------------------
    # 1. TRANSFORM SELF-TEST
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. INDEPENDENT FALLING-BASIS TRANSFORM SELF-TEST")
    print("=" * 78)

    (
        transform_ok,
        bad_d,
        bad_r,
        got,
        expected,
    ) = transform_self_test()

    print(f"transform_self_test = {transform_ok}")

    if not transform_ok:
        print(f"  failure degree={bad_d} r={bad_r}")
        print(f"  got={got}")
        print(f"  expected={expected}")

    print()

    # ------------------------------------------------------------------------
    # 2. BUILD Q_p(k)
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("2. EXACT TERMINAL-QUOTIENT EXTRACTION")
    print("=" * 78)

    quotients = {}
    extraction_ok = True

    for name, cfg in CHANNELS.items():
        K = cfg["K"]
        c = cfg["c"]

        quotients[name] = {}

        print(name)

        for p, values in cfg["data"].items():
            s = support_s(p, c)

            C = interpolate(values)
            D = terminal(K, s)

            try:
                Q = div_exact(C, D)
                exact = True
            except Exception:
                Q = [Fraction(0)]
                exact = False
                extraction_ok = False

            quotients[name][p] = Q

            print(
                f"  p={p}: s={s} "
                f"deg(C)={degree(C)} "
                f"deg(Q)={degree(Q)} "
                f"exact={exact}"
            )

        print()

    # ------------------------------------------------------------------------
    # 3. CORRECT SECOND FALLING-BASIS RECONSTRUCTION
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("3. CORRECT SECOND FALLING-BASIS RECONSTRUCTION")
    print("=" * 78)

    reconstruction_ok = True
    falling_rows = {}

    for name, cfg in CHANNELS.items():
        print(name)

        falling_rows[name] = {}

        for p, Q in quotients[name].items():
            coeffs = to_falling_basis(Q)
            rebuilt = from_falling_basis(coeffs)

            exact = rebuilt == trim(Q)

            if not exact:
                reconstruction_ok = False

            falling_rows[name][p] = coeffs

            print(f"  p={p}:")
            print(f"    q = {coeffs}")
            print(f"    reconstruction={exact}")

        print()

    # ------------------------------------------------------------------------
    # 4. POINTWISE INDEPENDENT CHECK
    #
    # This confirms that the coefficient vector represents Q at every
    # integer point, independently of polynomial coefficient comparison.
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("4. POINTWISE FALLING-BASIS CHECK")
    print("=" * 78)

    pointwise_ok = True

    for name, cfg in CHANNELS.items():
        K = cfg["K"]

        print(name)

        for p, coeffs in falling_rows[name].items():
            ok = True

            Q = quotients[name][p]

            for k in range(K + 1):
                lhs = eval_poly(Q, k)

                rhs = Fraction(0)

                for r, q in enumerate(coeffs):
                    rhs += q * Fraction(
                        factorial(k),
                        factorial(k - r),
                    ) if r <= k else Fraction(0)

                if lhs != rhs:
                    ok = False
                    break

            if not ok:
                pointwise_ok = False

            print(f"  p={p}: exact={ok}")

        print()

    # ------------------------------------------------------------------------
    # 5. SECOND FALLING-BASIS SIGNATURES
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("5. SECOND FALLING-BASIS PRIMITIVE SIGNATURES")
    print("=" * 78)

    for name, rows in falling_rows.items():
        print(name)

        for p, coeffs in rows.items():
            print(
                f"  p={p}: "
                f"{primitive_signature(coeffs)}"
            )

        print()

    # ------------------------------------------------------------------------
    # 6. SECOND-LAYER MATRIX RANK
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("6. SECOND-LAYER MATRIX RANK")
    print("=" * 78)

    ranks = {}

    for name, rows in falling_rows.items():
        matrix = list(rows.values())

        width = max(len(row) for row in matrix)

        matrix = [
            row + [Fraction(0)] * (width - len(row))
            for row in matrix
        ]

        ranks[name] = rank(matrix)

        print(
            f"{name}: shape=({len(matrix)}, {width}) "
            f"rank={ranks[name]}"
        )

    print()

    # ------------------------------------------------------------------------
    # 7. SECOND-LAYER TERMINAL / LOWER FACTORS
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("7. SECOND-LAYER EXACT FACTOR SEARCH")
    print("=" * 78)

    for name, rows in quotients.items():
        print(name)

        for p, Q in rows.items():
            hits = exact_factor_hits(Q)

            print(
                f"  p={p}: "
                f"{hits}"
            )

        print()

    # ------------------------------------------------------------------------
    # 8. SECOND-LAYER COEFFICIENT SUPPORT
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("8. SECOND-LAYER SUPPORT PROFILE")
    print("=" * 78)

    support_ok = True

    for name, rows in falling_rows.items():
        print(name)

        for p, coeffs in rows.items():
            nonzero = [
                r for r, x in enumerate(coeffs)
                if x != 0
            ]

            if nonzero:
                first = min(nonzero)
                last = max(nonzero)
            else:
                first = None
                last = None

            # The expected last index is simply deg(Q), because
            # k_(r) has degree r and the basis is triangular.
            expected_last = degree(quotients[name][p])

            exact = last == expected_last

            if not exact:
                support_ok = False

            print(
                f"  p={p}: "
                f"first={first} last={last} "
                f"degree={expected_last} exact={exact}"
            )

        print()

    # ------------------------------------------------------------------------
    # 9. CROSS-CHANNEL CORRECTED BASIS COMPARISON
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("9. CROSS-CHANNEL SECOND-BASIS COMPARISON")
    print("=" * 78)

    common_p = [
        0, 1, 2, 3, 4, 5, 6, 7
    ]

    for p in common_p:
        if (
            p not in falling_rows["A-even"]
            and p not in falling_rows["B-even"]
        ):
            continue

        pairs = []

        if p in falling_rows["A-even"] and p in falling_rows["B-even"]:
            pairs.append(("A-even", "B-even"))

        if p in falling_rows["A-odd"] and p in falling_rows["B-odd"]:
            pairs.append(("A-odd", "B-odd"))

        for left, right in pairs:
            A = falling_rows[left][p]
            B = falling_rows[right][p]

            width = max(len(A), len(B))
            A = A + [Fraction(0)] * (width - len(A))
            B = B + [Fraction(0)] * (width - len(B))

            ratios = []

            for r in range(width):
                if B[r] != 0:
                    ratios.append(
                        (r, A[r] / B[r])
                    )
                elif A[r] != 0:
                    ratios.append(
                        (r, None)
                    )

            print(
                f"p={p}: {left} / {right}"
            )
            print(f"  ratios={ratios}")

        if pairs:
            print()

    # ------------------------------------------------------------------------
    # 10. NEW QUESTION:
    #
    # Does the CORRECT residual falling-basis matrix admit a triangular
    # relation with the original coefficient index p?
    #
    # Since Q_p is already after removal of the terminal factor, this tests
    # whether a relation of the rough form
    #
    #     q[p,r] = 0 for r < h(p)
    #
    # exists for any simple h(p).
    #
    # We only report the exact observed support; we do NOT extrapolate.
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("10. INDEX-SUPPORT AUDIT IN THE CORRECTED BASIS")
    print("=" * 78)

    for name, rows in falling_rows.items():
        print(name)

        for p, coeffs in rows.items():
            nz = [
                r for r, x in enumerate(coeffs)
                if x != 0
            ]

            print(
                f"  p={p}: "
                f"support={nz}"
            )

        print()

    # ------------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------------

    failures = 0

    checks = {
        "transform_self_test": transform_ok,
        "terminal_extraction": extraction_ok,
        "falling_reconstruction": reconstruction_ok,
        "pointwise_falling_check": pointwise_ok,
        "second_basis_support": support_ok,
    }

    for ok in checks.values():
        if not ok:
            failures += 1

    print("=" * 78)
    print("11. STRUCTURAL INTERPRETATION")
    print("=" * 78)
    print()
    print("Experiment 114's apparent failure was computational, not")
    print("structural: its falling-basis transform used floor division")
    print("on exact Fractions.")
    print()
    print("This experiment repairs that operation and independently")
    print("tests the basis identity:")
    print()
    print("    Q_p(k) = sum_r q[p,r] k_(r).")
    print()
    print("The key distinction is now clean:")
    print()
    print("    FIRST FACTOR:")
    print("        C_p(k) = (K-k)_s Q_p(k)")
    print()
    print("    SECOND BASIS:")
    print("        Q_p(k) = sum_r q[p,r] k_(r)")
    print()
    print("A genuinely new structural discovery would require the corrected")
    print("q[p,r] matrix to exhibit additional exact triangularity, factors,")
    print("low rank, or a simple cross-channel relation.")
    print()
    print("No such conclusion is assumed in advance.")
    print("Everything is exact over QQ.")
    print("No SymPy.")
    print("No floating point.")
    print("No extrapolation.")
    print()

    print("=" * 78)
    print("12. FINAL EXACTNESS")
    print("=" * 78)

    print(f"  transform_self_test = {transform_ok}")
    print(f"  terminal_extraction = {extraction_ok}")
    print(f"  falling_reconstruction = {reconstruction_ok}")
    print(f"  pointwise_falling_check = {pointwise_ok}")
    print(f"  second_basis_support = {support_ok}")
    print(f"  failures = {failures}")
    print(
        f"  ALL BASIC CHECKS PASS = {failures == 0}"
    )

    print()
    print("EXPERIMENT 115 COMPLETE")


if __name__ == "__main__":
    main()

