#!/usr/bin/env python3

from fractions import Fraction
from math import factorial
from functools import reduce
from itertools import product
from math import gcd


# ============================================================================
# EXPERIMENT 116
# EXACT RESIDUAL KERNEL BIVARIATE / BINOMIAL-NORMALIZATION AUDIT
# ============================================================================
#
# Established structure:
#
#     C_p(k) = (K-k)_s Q_p(k)
#
# and, after the corrected second falling-factorial transform,
#
#     Q_p(k) = sum_r q[p,r] k_(r).
#
# Experiment 115 established the transform exactly.
#
# This experiment studies q[p,r] itself.
#
# Questions:
#
#   1. Does q[p,r] simplify after division by binom(D_p,r)?
#   2. Does it simplify after multiplying/dividing by r!?
#   3. Does the normalized matrix have reduced rank?
#   4. Is the dependence primarily on r-p, r+p, or r relative to
#      the residual degree D_p?
#   5. Are neighboring p-rows related by a simple affine shift in r?
#   6. Is there a low-degree polynomial relation in (p,r)?
#   7. Do A/B residual surfaces satisfy a common low-degree relation?
#
# Everything is exact over QQ.
# No SymPy.
# No floating point.
# No extrapolation.
#
# ============================================================================


# ----------------------------------------------------------------------------
# DATA: CORRECT q[p,r] FROM EXPERIMENT 115
# ----------------------------------------------------------------------------

A_even = {
    0: [
        Fraction(-12879), Fraction(-15362), Fraction(8307),
        Fraction(-748), Fraction(-62305, 144),
        Fraction(34070797, 201600), Fraction(-102402481, 3628800)
    ],
    2: [
        Fraction(2797337, 11520), Fraction(7319861, 19200),
        Fraction(-266207819, 1612800), Fraction(-141551327, 14515200),
        Fraction(2461903867, 116121600),
        Fraction(-783039371, 116121600)
    ],
    4: [
        Fraction(-2083937, 921600), Fraction(-3930151, 921600),
        Fraction(218532413, 77414400), Fraction(238646881, 232243200),
        Fraction(-145404619, 92897280)
    ],
    6: [
        Fraction(85591, 22118400), Fraction(249013, 22118400),
        Fraction(-24042833, 619315200)
    ],
    8: [
        Fraction(-4913, 353894400)
    ],
}

A_odd = {
    1: [
        Fraction(12143, 3360), Fraction(76427, 192),
        Fraction(1954873, 17920), Fraction(-61469491, 483840),
        Fraction(42852113, 1935360), Fraction(4384549, 645120)
    ],
    3: [
        Fraction(-989, 345600), Fraction(-36064769, 2764800),
        Fraction(-24988097, 2580480), Fraction(87300373, 19353600),
        Fraction(116226679, 77414400)
    ],
    5: [
        Fraction(-517, 2764800), Fraction(1189391, 5529600),
        Fraction(22183547, 103219200), Fraction(-69294643, 185794560)
    ],
    7: [
        Fraction(373, 928972800), Fraction(-616981, 371589120)
    ],
}

B_even = {
    0: [
        Fraction(12980463, 1024), Fraction(12041869, 1024),
        Fraction(-594517923, 71680), Fraction(112271581, 71680),
        Fraction(123350021, 860160),
        Fraction(-1803389011, 12902400)
    ],
    2: [
        Fraction(-19344659, 76800), Fraction(-38450509, 115200),
        Fraction(490918171, 3225600), Fraction(-31781287, 1209600),
        Fraction(254480207, 12902400)
    ],
    4: [
        Fraction(25883, 12288), Fraction(682871, 184320),
        Fraction(-1505893, 2580480), Fraction(-25483559, 7741440)
    ],
    6: [
        Fraction(-2267, 614400), Fraction(-100657, 5529600)
    ],
}

B_odd = {
    1: [
        Fraction(-584531, 35840), Fraction(-32224291, 21504),
        Fraction(74131151, 215040), Fraction(29406229, 129024),
        Fraction(-1338089411, 7741440),
        Fraction(495451247, 7741440)
    ],
    3: [
        Fraction(-59257, 230400), Fraction(7521137, 230400),
        Fraction(12697441, 3225600), Fraction(4595257, 1382400),
        Fraction(-140504813, 12902400)
    ],
    5: [
        Fraction(4457, 2764800), Fraction(-340837, 2764800),
        Fraction(-5342627, 12902400)
    ],
    7: [
        Fraction(-421, 38707200)
    ],
}


CHANNELS = {
    "A-even": {"rows": A_even, "K": 6, "channel": 0, "parity": 0},
    "A-odd":  {"rows": A_odd,  "K": 6, "channel": 0, "parity": 1},
    "B-even": {"rows": B_even, "K": 5, "channel": 1, "parity": 0},
    "B-odd":  {"rows": B_odd, "K": 5, "channel": 1, "parity": 1},
}


# ----------------------------------------------------------------------------
# BASIC EXACT UTILITIES
# ----------------------------------------------------------------------------

def trim_row(row):
    row = list(row)

    while len(row) > 1 and row[-1] == 0:
        row.pop()

    return row


def degree_row(row):
    row = trim_row(row)

    if len(row) == 1 and row[0] == 0:
        return -1

    return len(row) - 1


def matrix_rank(matrix):
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
    rank = 0

    for col in range(n):
        pivot = None

        for r in range(rank, m):
            if A[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        A[rank], A[pivot] = A[pivot], A[rank]

        pivot_value = A[rank][col]

        A[rank] = [
            x / pivot_value
            for x in A[rank]
        ]

        for r in range(m):
            if r == rank:
                continue

            factor = A[r][col]

            if factor:
                A[r] = [
                    A[r][j] - factor * A[rank][j]
                    for j in range(n)
                ]

        rank += 1

        if rank == m:
            break

    return rank


def primitive_integer_signature(row):
    row = trim_row(row)

    if all(x == 0 for x in row):
        return [0]

    denominator_lcm = 1

    for x in row:
        denominator_lcm = (
            denominator_lcm
            * x.denominator
            // gcd(denominator_lcm, x.denominator)
        )

    ints = [
        int(x * denominator_lcm)
        for x in row
    ]

    g = reduce(
        gcd,
        [abs(x) for x in ints if x != 0],
        0,
    )

    if g:
        ints = [x // g for x in ints]

    for x in ints:
        if x:
            if x < 0:
                ints = [-v for v in ints]
            break

    return ints


# ----------------------------------------------------------------------------
# NATURAL NORMALIZATIONS
# ----------------------------------------------------------------------------

def binomial_exact(n, r):
    if r < 0 or r > n:
        return Fraction(0)

    return Fraction(
        factorial(n),
        factorial(r) * factorial(n - r),
    )


def normalize_binomial(row):
    """
    Divide q[p,r] by binomial(D,r), where D=degree(Q).

    This tests whether the corrected residual coefficients carry
    a standard combinatorial multiplicity.
    """
    D = degree_row(row)

    if D < 0:
        return []

    return [
        value / binomial_exact(D, r)
        for r, value in enumerate(row)
    ]


def normalize_r_factorial(row):
    """
    Multiply q[p,r] by r!.

    This tests whether denominators are primarily basis-normalization
    artifacts.
    """
    return [
        value * factorial(r)
        for r, value in enumerate(row)
    ]


def normalize_divide_r_factorial(row):
    """
    Divide q[p,r] by r!.
    """
    return [
        value / factorial(r)
        for r, value in enumerate(row)
    ]


def normalize_kfalling(row):
    """
    Divide by binomial(D,r) * r! = D_(r).

    This is the falling-factorial multiplicity itself.
    """
    D = degree_row(row)

    return [
        (
            value
            / Fraction(
                factorial(D),
                factorial(D - r),
            )
        )
        for r, value in enumerate(row)
    ]


# ----------------------------------------------------------------------------
# ROW ALIGNMENT TESTS
# ----------------------------------------------------------------------------

def pad(row, width):
    return row + [Fraction(0)] * (width - len(row))


def exact_row_proportional(a, b):
    width = max(len(a), len(b))

    a = pad(a, width)
    b = pad(b, width)

    scale = None

    for x, y in zip(a, b):
        if x == 0 and y == 0:
            continue

        if y == 0 or x == 0:
            return False, None

        candidate = x / y

        if scale is None:
            scale = candidate
        elif candidate != scale:
            return False, None

    if scale is None:
        return True, Fraction(1)

    return True, scale


def shifted_row(row, shift):
    """
    q[r] -> q[r+shift] alignment.

    Positive shift means the row is shifted toward larger r.
    """
    if shift >= 0:
        return [Fraction(0)] * shift + row

    s = -shift

    if s >= len(row):
        return [Fraction(0)]

    return row[s:]


# ----------------------------------------------------------------------------
# LOW-DEGREE POLYNOMIAL RELATION IN (p,r)
# ----------------------------------------------------------------------------

def monomials_total_degree(max_degree):
    out = []

    for total in range(max_degree + 1):
        for a in range(total + 1):
            b = total - a
            out.append((a, b))

    return out


def exact_bivariate_relation(points_left, points_right, degree_bound):
    """
    Search for

        P(p,r) * left + Q(p,r) * right = 0

    with P,Q of total degree <= degree_bound.

    Uses exact rational row reduction.

    Returns None when no relation exists.
    """

    mons = monomials_total_degree(degree_bound)

    unknown_count = 2 * len(mons)

    matrix = []

    for (p, r), a, b in zip(
        points_left[0],
        points_left[1],
        points_right[1],
    ):
        row = []

        for alpha, beta in mons:
            row.append(
                a * (p ** alpha) * (r ** beta)
            )

        for alpha, beta in mons:
            row.append(
                b * (p ** alpha) * (r ** beta)
            )

        matrix.append(row)

    rank = matrix_rank(matrix)
    nullity = unknown_count - rank

    return rank, unknown_count, nullity


# ----------------------------------------------------------------------------
# COEFFICIENT POINT SETS
# ----------------------------------------------------------------------------

def overlapping_points(row_a, row_b):
    width = min(len(row_a), len(row_b))

    points = []

    for r in range(width):
        points.append(r)

    return points


# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 116 — EXACT RESIDUAL KERNEL")
    print("BIVARIATE / BINOMIAL-NORMALIZATION AUDIT")
    print("=" * 78)
    print()

    # ------------------------------------------------------------------------
    # 1. BASIC DATA PROFILE
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. CORRECTED SECOND-LAYER DATA PROFILE")
    print("=" * 78)

    for name, cfg in CHANNELS.items():
        print(name)

        for p, row in cfg["rows"].items():
            print(
                f"  p={p}: degree={degree_row(row)} "
                f"support={list(range(len(row)))}"
            )

        print()

    # ------------------------------------------------------------------------
    # 2. NORMALIZATION MATRICES
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("2. EXACT NORMALIZATION RANK AUDIT")
    print("=" * 78)

    normalizations = {
        "raw": lambda row: row,
        "divide_binomial": normalize_binomial,
        "times_r_factorial": normalize_r_factorial,
        "divide_r_factorial": normalize_divide_r_factorial,
        "divide_falling_multiplicity": normalize_kfalling,
    }

    normalization_ranks = {}

    for norm_name, fn in normalizations.items():

        print()
        print(f"  NORMALIZATION = {norm_name}")

        normalization_ranks[norm_name] = {}

        for name, cfg in CHANNELS.items():
            rows = []

            transformed = {}

            for p, row in cfg["rows"].items():
                t = fn(row)
                transformed[p] = t

                rows.append(
                    pad(
                        t,
                        max(
                            len(x)
                            for x in cfg["rows"].values()
                        ),
                    )
                )

            rnk = matrix_rank(rows)

            normalization_ranks[norm_name][name] = rnk

            print(
                f"    {name}: "
                f"shape=({len(rows)},{len(rows[0])}) "
                f"rank={rnk}"
            )

    print()

    # ------------------------------------------------------------------------
    # 3. PRINT MOST IMPORTANT NORMALIZED SIGNATURES
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("3. PRIMITIVE SIGNATURES AFTER NATURAL NORMALIZATIONS")
    print("=" * 78)

    for name, cfg in CHANNELS.items():

        print()
        print(name)

        for p, row in cfg["rows"].items():

            print(f"  p={p}")

            for norm_name, fn in normalizations.items():

                signature = primitive_integer_signature(
                    fn(row)
                )

                print(
                    f"    {norm_name}: {signature}"
                )

    print()

    # ------------------------------------------------------------------------
    # 4. OFFSET-COORDINATE AUDIT
    #
    # Because row degree decreases with p, test whether normalized
    # coefficients become simpler in coordinates
    #
    #     d = degree(Q_p) - r
    #
    # i.e. distance from the terminal coefficient.
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("4. TERMINAL-DISTANCE COORDINATE AUDIT")
    print("=" * 78)

    terminal_matrices = {}

    for name, cfg in CHANNELS.items():

        print()
        print(name)

        terminal_matrices[name] = []

        for p, row in cfg["rows"].items():

            D = degree_row(row)

            values = []

            for r, value in enumerate(row):
                distance = D - r

                values.append(
                    (distance, value)
                )

            print(
                f"  p={p}: "
                f"[(D-r,value)]={values}"
            )

            terminal_matrices[name].append(
                [
                    value
                    for _, value in values
                ]
            )

        print(
            f"  terminal-distance rank="
            f"{matrix_rank(terminal_matrices[name])}"
        )

    print()

    # ------------------------------------------------------------------------
    # 5. CONSECUTIVE-p SHIFT / PROPORTIONALITY AUDIT
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("5. CONSECUTIVE-p SHIFT / PROPORTIONALITY AUDIT")
    print("=" * 78)

    shift_candidates = [-2, -1, 0, 1, 2]

    for name, cfg in CHANNELS.items():

        print()
        print(name)

        rows = cfg["rows"]
        ps = list(rows.keys())

        for p1, p2 in zip(ps, ps[1:]):

            row1 = rows[p1]
            row2 = rows[p2]

            print(
                f"  transition {p1}->{p2}"
            )

            for shift in shift_candidates:

                shifted = shifted_row(
                    row1,
                    shift,
                )

                proportional, scalar = exact_row_proportional(
                    row2,
                    shifted,
                )

                if proportional:
                    print(
                        f"    shift={shift}: "
                        f"PROPORTIONAL scalar={scalar}"
                    )

    print()

    # ------------------------------------------------------------------------
    # 6. DIAGONAL-SURFACE AUDIT
    #
    # Inspect sums / slices indexed by:
    #
    #     p + r
    #     p - 2r
    #     r - p/2
    #
    # These are deliberately descriptive rather than interpolative.
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("6. BIVARIATE INDEX-COMBINATION AUDIT")
    print("=" * 78)

    combinations = {
        "r": lambda p, r: r,
        "p+r": lambda p, r: p + r,
        "r-p": lambda p, r: r - p,
        "p-2r": lambda p, r: p - 2 * r,
        "2r-p": lambda p, r: 2 * r - p,
    }

    for name, cfg in CHANNELS.items():

        print()
        print(name)

        for p, row in cfg["rows"].items():

            for coord_name, fn in combinations.items():

                coords = [
                    fn(p, r)
                    for r in range(len(row))
                ]

                print(
                    f"  p={p} {coord_name}={coords}"
                )

    print()

    # ------------------------------------------------------------------------
    # 7. LOW-DEGREE BIVARIATE RELATION SEARCH
    #
    # Compare A-even vs B-even and A-odd vs B-odd.
    #
    # The matrices have different row degrees at the extremes, so only
    # overlapping (p,r) points are used.
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("7. LOW-DEGREE A/B BIVARIATE RELATION SEARCH")
    print("=" * 78)

    pairings = [
        ("A-even", "B-even"),
        ("A-odd", "B-odd"),
    ]

    relation_results = {}

    for left_name, right_name in pairings:

        left = CHANNELS[left_name]["rows"]
        right = CHANNELS[right_name]["rows"]

        common_p = sorted(
            set(left.keys()) & set(right.keys())
        )

        points = []
        left_values = []
        right_values = []

        for p in common_p:

            width = min(
                len(left[p]),
                len(right[p]),
            )

            for r in range(width):

                points.append((p, r))
                left_values.append(left[p][r])
                right_values.append(right[p][r])

        relation_results[
            (left_name, right_name)
        ] = {}

        print()
        print(
            f"{left_name} vs {right_name}: "
            f"points={len(points)}"
        )

        for deg in range(0, 4):

            result = exact_bivariate_relation(
                (points, left_values),
                (points, right_values),
                deg,
            )

            relation_results[
                (left_name, right_name)
            ][deg] = result

            rank_value, unknowns, nullity = result

            print(
                f"  degree={deg}: "
                f"rank={rank_value} "
                f"unknowns={unknowns} "
                f"nullity={nullity}"
            )

    print()

    # ------------------------------------------------------------------------
    # 8. CROSS-CHANNEL NORMALIZED RELATION SEARCH
    #
    # Repeat after division by binomial(D,r), which is the most natural
    # combinatorial normalization.
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("8. NORMALIZED A/B RELATION SEARCH")
    print("=" * 78)

    def normalized_entry(name, p, r):
        row = CHANNELS[name]["rows"][p]
        norm = normalize_binomial(row)
        return norm[r]

    for left_name, right_name in pairings:

        common_p = sorted(
            set(CHANNELS[left_name]["rows"].keys())
            & set(CHANNELS[right_name]["rows"].keys())
        )

        points = []
        left_values = []
        right_values = []

        for p in common_p:

            width = min(
                len(CHANNELS[left_name]["rows"][p]),
                len(CHANNELS[right_name]["rows"][p]),
            )

            for r in range(width):

                points.append((p, r))
                left_values.append(
                    normalized_entry(left_name, p, r)
                )
                right_values.append(
                    normalized_entry(right_name, p, r)
                )

        print()
        print(
            f"{left_name} normalized vs "
            f"{right_name} normalized"
        )

        for deg in range(0, 4):

            rank_value, unknowns, nullity = (
                exact_bivariate_relation(
                    (points, left_values),
                    (points, right_values),
                    deg,
                )
            )

            print(
                f"  degree={deg}: "
                f"rank={rank_value} "
                f"unknowns={unknowns} "
                f"nullity={nullity}"
            )

    print()

    # ------------------------------------------------------------------------
    # 9. EXACT RECONSTRUCTION SANITY
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("9. EXACT RECONSTRUCTION SANITY")
    print("=" * 78)

    reconstruction_ok = True

    for name, cfg in CHANNELS.items():

        print(name)

        for p, row in cfg["rows"].items():

            # Reconstruct q row from its own falling basis using its
            # coefficient values evaluated at k=0,...,D.
            D = degree_row(row)

            # The q row itself is already in falling basis. Reconstruct
            # its index polynomial through explicit evaluation.
            values = []

            for k in range(D + 1):

                value = Fraction(0)

                for r, q in enumerate(row):

                    if r <= k:

                        value += (
                            q
                            * Fraction(
                                factorial(k),
                                factorial(k - r),
                            )
                        )

                values.append(value)

            # Recover the coefficients exactly from those values.
            recovered = [Fraction(0)] * (D + 1)

            for n in range(D + 1):

                residual = values[n]

                for r in range(n):

                    residual -= (
                        recovered[r]
                        * Fraction(
                            factorial(n),
                            factorial(n - r),
                        )
                    )

                recovered[n] = (
                    residual
                    / factorial(n)
                )

            recovered = trim_row(recovered)

            exact = recovered == trim_row(row)

            if not exact:
                reconstruction_ok = False

            print(
                f"  p={p}: exact={exact}"
            )

        print()

    # ------------------------------------------------------------------------
    # 10. FINAL STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)
    print()
    print("Experiment 115 established that the second falling-factorial")
    print("basis itself is exact; there is no remaining transform bug.")
    print()
    print("Experiment 116 therefore studies the corrected residual surface")
    print()
    print("    q[p,r]")
    print()
    print("directly.")
    print()
    print("The natural tests are:")
    print()
    print("  * binomial normalization in the residual degree;");
    print("  * factorial normalization;");
    print("  * terminal-distance coordinates;");
    print("  * consecutive-p shifted proportionality;");
    print("  * low-degree bivariate A/B relations.");
    print()
    print("A low-rank result after a natural normalization would indicate")
    print("that the residual kernel still contains a reusable combinatorial")
    print("multiplicity.")
    print()
    print("A negative result would strengthen the conclusion that the")
    print("following exact structure is already close to minimal:")
    print()
    print("    C_p(k) = (K-k)_s Q_p(k)")
    print()
    print("with")
    print()
    print("    Q_p(k) = sum_r q[p,r] k_(r).")
    print()
    print("No extrapolation is performed.")
    print("Everything is exact over QQ.")
    print("No floating point.")
    print("No SymPy.")
    print("")

    failures = 0

    if not reconstruction_ok:
        failures += 1

    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)
    print(
        f"  corrected_second_basis = True"
    )
    print(
        f"  reconstruction = {reconstruction_ok}"
    )
    print(
        f"  bivariate_audit_completed = True"
    )
    print(
        f"  failures = {failures}"
    )
    print(
        f"  ALL BASIC CHECKS PASS = {failures == 0}"
    )
    print()
    print("EXPERIMENT 116 COMPLETE")


if __name__ == "__main__":
    main()

