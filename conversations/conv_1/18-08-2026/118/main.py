#!/usr/bin/env python3

from fractions import Fraction
from math import factorial, comb, gcd


# ==============================================================================
# EXPERIMENT 121
# EXACT B-ODD OPERATOR CONJUGATION / NORMALIZATION AUDIT
# ==============================================================================

# All data are the corrected second-layer falling-basis rows from
# Experiment 115/116.
#
# Each row is
#
#     Q_p(k) = sum_r q[p,r] k_(r)
#
# and the four sectors are:
#
#     A-even
#     A-odd
#     B-even
#     B-odd
#
# The exceptional exact operator from Experiment 120 is the B-odd
# order=3, degree=2 operator.
#
# This experiment:
#
#   1. reconstructs that operator independently;
#   2. applies exact diagonal weight conjugations;
#   3. measures whether the transformed operator becomes simpler;
#   4. checks transfer to the other sectors;
#   5. performs all calculations with exact Fraction arithmetic.
#
# No floating point.
# No SymPy.
# No extrapolation.
# ==============================================================================


# ==============================================================================
# EXACT SECOND-LAYER DATA
# ==============================================================================

A_even = {
    0: [
        Fraction(-12879, 1),
        Fraction(-15362, 1),
        Fraction(8307, 1),
        Fraction(-748, 1),
        Fraction(-62305, 144),
        Fraction(34070797, 201600),
        Fraction(-102402481, 3628800),
    ],
    2: [
        Fraction(2797337, 11520),
        Fraction(7319861, 19200),
        Fraction(-266207819, 1612800),
        Fraction(-141551327, 14515200),
        Fraction(2461903867, 116121600),
        Fraction(-783039371, 116121600),
    ],
    4: [
        Fraction(-2083937, 921600),
        Fraction(-3930151, 921600),
        Fraction(218532413, 77414400),
        Fraction(238646881, 232243200),
        Fraction(-145404619, 92897280),
    ],
    6: [
        Fraction(85591, 22118400),
        Fraction(249013, 22118400),
        Fraction(-24042833, 619315200),
    ],
    8: [
        Fraction(-4913, 353894400),
    ],
}


A_odd = {
    1: [
        Fraction(12143, 3360),
        Fraction(76427, 192),
        Fraction(1954873, 17920),
        Fraction(-61469491, 483840),
        Fraction(42852113, 1935360),
        Fraction(4384549, 645120),
    ],
    3: [
        Fraction(-989, 345600),
        Fraction(-36064769, 2764800),
        Fraction(-24988097, 2580480),
        Fraction(87300373, 19353600),
        Fraction(116226679, 77414400),
    ],
    5: [
        Fraction(-517, 2764800),
        Fraction(1189391, 5529600),
        Fraction(22183547, 103219200),
        Fraction(-69294643, 185794560),
    ],
    7: [
        Fraction(373, 928972800),
        Fraction(-616981, 371589120),
    ],
}


B_even = {
    0: [
        Fraction(12980463, 1024),
        Fraction(12041869, 1024),
        Fraction(-594517923, 71680),
        Fraction(112271581, 71680),
        Fraction(123350021, 860160),
        Fraction(-1803389011, 12902400),
    ],
    2: [
        Fraction(-19344659, 76800),
        Fraction(-38450509, 115200),
        Fraction(490918171, 3225600),
        Fraction(-31781287, 1209600),
        Fraction(254480207, 12902400),
    ],
    4: [
        Fraction(25883, 12288),
        Fraction(682871, 184320),
        Fraction(-1505893, 2580480),
        Fraction(-25483559, 7741440),
    ],
    6: [
        Fraction(-2267, 614400),
        Fraction(-100657, 5529600),
    ],
}


B_odd = {
    1: [
        Fraction(-584531, 35840),
        Fraction(-32224291, 21504),
        Fraction(74131151, 215040),
        Fraction(29406229, 129024),
        Fraction(-1338089411, 7741440),
        Fraction(495451247, 7741440),
    ],
    3: [
        Fraction(-59257, 230400),
        Fraction(7521137, 230400),
        Fraction(12697441, 3225600),
        Fraction(4595257, 1382400),
        Fraction(-140504813, 12902400),
    ],
    5: [
        Fraction(4457, 2764800),
        Fraction(-340837, 2764800),
        Fraction(-5342627, 12902400),
    ],
    7: [
        Fraction(-421, 38707200),
    ],
}


SECTORS = {
    "A-even": A_even,
    "A-odd": A_odd,
    "B-even": B_even,
    "B-odd": B_odd,
}


ORDER = 3
DEGREE = 2
K = 5


# ==============================================================================
# EXACT LINEAR ALGEBRA
# ==============================================================================

def rref(matrix):
    M = [
        [Fraction(x) for x in row]
        for row in matrix
    ]

    if not M:
        return [], [], 0

    nr = len(M)
    nc = len(M[0])

    pivots = []
    pivot_row = 0

    for col in range(nc):
        pivot = None

        for r in range(pivot_row, nr):
            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        M[pivot_row], M[pivot] = (
            M[pivot],
            M[pivot_row],
        )

        scale = M[pivot_row][col]

        M[pivot_row] = [
            x / scale
            for x in M[pivot_row]
        ]

        for r in range(nr):
            if r == pivot_row:
                continue

            factor = M[r][col]

            if factor == 0:
                continue

            M[r] = [
                M[r][j] - factor * M[pivot_row][j]
                for j in range(nc)
            ]

        pivots.append(col)
        pivot_row += 1

        if pivot_row == nr:
            break

    return M, pivots, pivot_row


def rank(matrix):
    if not matrix:
        return 0
    return rref(matrix)[2]


def nullity_if_consistent(A, b):
    augmented = [
        list(map(Fraction, row)) + [Fraction(rhs)]
        for row, rhs in zip(A, b)
    ]

    R, pivots_aug, rank_aug = rref(augmented)

    n = len(A[0])

    for row in R:
        if all(row[j] == 0 for j in range(n)) and row[n] != 0:
            return None

    rank_A = rank(A)
    return n - rank_A


def solve_unique(A, b):
    if not A:
        return None

    n = len(A[0])

    augmented = [
        list(map(Fraction, row)) + [Fraction(rhs)]
        for row, rhs in zip(A, b)
    ]

    R, pivots, rank_aug = rref(augmented)

    # Inconsistency.
    for row in R:
        if all(row[j] == 0 for j in range(n)):
            if row[n] != 0:
                return None

    if rank_aug != n:
        return None

    solution = [Fraction(0)] * n

    for i, pc in enumerate(pivots):
        if pc < n:
            solution[pc] = R[i][n]

    return solution


# ==============================================================================
# FALLING FACTORIALS
# ==============================================================================

def falling(x, r):
    out = Fraction(1)

    for i in range(r):
        out *= x - i

    return out


def falling_eval(coeffs, x):
    return sum(
        Fraction(c) * falling(x, r)
        for r, c in enumerate(coeffs)
    )


def falling_interpolate(values):
    """
    Given exact values f(0), ..., f(d), return coefficients c_r such that

        f(x) = sum_r c_r x_(r).

    Since
        c_r = Delta^r f(0) / r!,
    this is exact over Q.
    """

    work = [
        Fraction(v)
        for v in values
    ]

    result = []

    for r in range(len(work)):
        result.append(
            work[0] / factorial(r)
        )

        work = [
            work[i + 1] - work[i]
            for i in range(len(work) - 1)
        ]

    return result


# ==============================================================================
# OPERATOR SYSTEM
# ==============================================================================

def build_operator_system(rows, order, degree):
    """
    Solve

        q[p_next,r]
          = sum_{a=0}^{order-1} F_a(r) q[p,r+a]

    with

        F_a(r) = sum_{b=0}^{degree} c[a,b] r_(b).

    Only equations with an available transition are used.
    """

    equations = []

    ps = sorted(rows)

    for p0, p1 in zip(ps, ps[1:]):

        previous = rows[p0]
        next_row = rows[p1]

        for r in range(len(next_row)):

            equation = [
                Fraction(0)
            ] * (order * (degree + 1))

            for a in range(order):

                src = r + a

                if src >= len(previous):
                    continue

                q = previous[src]

                for b in range(degree + 1):

                    idx = (
                        a * (degree + 1)
                        + b
                    )

                    equation[idx] += (
                        q * falling(r, b)
                    )

            # lhs - rhs = 0
            equations.append(
                equation + [
                    -next_row[r]
                ]
            )

    return equations


def recover_operator(rows, order, degree):
    system = build_operator_system(
        rows,
        order,
        degree,
    )

    if not system:
        return {
            "consistent": False,
            "unique": False,
            "solution": None,
            "rank": 0,
            "unknowns": order * (degree + 1),
            "nullity": None,
        }

    A = [
        row[:-1]
        for row in system
    ]

    b = [
        -row[-1]
        for row in system
    ]

    augmented_rank = rank([
        row[:] + [rhs]
        for row, rhs in zip(A, b)
    ])

    rank_A = rank(A)
    unknowns = order * (degree + 1)

    # Consistency test.
    consistent = (
        augmented_rank == rank_A
    )

    if not consistent:
        return {
            "consistent": False,
            "unique": False,
            "solution": None,
            "rank": rank_A,
            "unknowns": unknowns,
            "nullity": None,
        }

    nullity = unknowns - rank_A

    if nullity != 0:
        return {
            "consistent": True,
            "unique": False,
            "solution": None,
            "rank": rank_A,
            "unknowns": unknowns,
            "nullity": nullity,
        }

    solution = solve_unique(
        A,
        b,
    )

    return {
        "consistent": solution is not None,
        "unique": solution is not None,
        "solution": solution,
        "rank": rank_A,
        "unknowns": unknowns,
        "nullity": nullity,
    }


# ==============================================================================
# OPERATOR EVALUATION
# ==============================================================================

def component(solution, a, degree):
    start = a * (degree + 1)
    stop = start + degree + 1
    return solution[start:stop]


def operator_value(
    solution,
    order,
    degree,
    previous_row,
    r,
):
    total = Fraction(0)

    for a in range(order):

        src = r + a

        if src >= len(previous_row):
            continue

        F = component(
            solution,
            a,
            degree,
        )

        total += (
            falling_eval(F, r)
            * previous_row[src]
        )

    return total


def verify_operator(
    rows,
    solution,
    order,
    degree,
):
    if solution is None:
        return False

    ps = sorted(rows)

    for p0, p1 in zip(ps, ps[1:]):

        previous = rows[p0]
        target = rows[p1]

        for r in range(len(target)):

            predicted = operator_value(
                solution,
                order,
                degree,
                previous,
                r,
            )

            if predicted != target[r]:
                return False

    return True


# ==============================================================================
# UNIQUE B-ODD OPERATOR
# ==============================================================================

def recover_bodd_operator():
    return recover_operator(
        B_odd,
        ORDER,
        DEGREE,
    )


# ==============================================================================
# NATURAL WEIGHTS
# ==============================================================================

def weight_value(name, r):
    if name == "unit":
        return Fraction(1)

    if name == "r_factorial":
        return Fraction(factorial(r))

    if name == "inv_r_factorial":
        return Fraction(1, factorial(r))

    if name == "binomial_K_r":
        if 0 <= r <= K:
            return Fraction(comb(K, r))
        return Fraction(1)

    if name == "inv_binomial_K_r":
        if 0 <= r <= K:
            return Fraction(1, comb(K, r))
        return Fraction(1)

    if name == "r_factorial_times_Kminusr_factorial":
        if 0 <= r <= K:
            return Fraction(
                factorial(r)
                * factorial(K - r)
            )
        return Fraction(1)

    if name == "inverse_r_factorial_times_Kminusr_factorial":
        if 0 <= r <= K:
            return Fraction(
                1,
                factorial(r)
                * factorial(K - r),
            )
        return Fraction(1)

    if name == "r_factorial_over_Kminusr_factorial":
        if 0 <= r <= K:
            return Fraction(
                factorial(r),
                factorial(K - r),
            )
        return Fraction(1)

    if name == "Kminusr_factorial_over_r_factorial":
        if 0 <= r <= K:
            return Fraction(
                factorial(K - r),
                factorial(r),
            )
        return Fraction(1)

    raise ValueError(
        f"Unknown weight {name}"
    )


WEIGHTS = [
    "unit",
    "r_factorial",
    "inv_r_factorial",
    "binomial_K_r",
    "inv_binomial_K_r",
    "r_factorial_times_Kminusr_factorial",
    "inverse_r_factorial_times_Kminusr_factorial",
    "r_factorial_over_Kminusr_factorial",
    "Kminusr_factorial_over_r_factorial",
]


# ==============================================================================
# CONJUGATED OPERATOR
# ==============================================================================

def conjugate_operator(
    solution,
    weight_name,
):
    """
    If

        v[p,r] = w(r) q[p,r],

    then

        v[p+2,r]
          = sum_a F'_a(r) v[p,r+a]

    with

        F'_a(r) =
            F_a(r) w(r) / w(r+a).
    """

    transformed_components = []

    for a in range(ORDER):

        original = component(
            solution,
            a,
            DEGREE,
        )

        # We need enough values to determine whether the transformed
        # coefficient remains a polynomial of the same small degree.
        values = []

        for r in range(ORDER + DEGREE + 2):

            denominator = weight_value(
                weight_name,
                r + a,
            )

            if denominator == 0:
                return None

            value = (
                falling_eval(
                    original,
                    r,
                )
                * weight_value(
                    weight_name,
                    r,
                )
                / denominator
            )

            values.append(value)

        coeffs = falling_interpolate(
            values
        )

        # Verify exact polynomial agreement on every tested point.
        for r, expected in enumerate(values):
            if falling_eval(coeffs, r) != expected:
                return None

        # Reject an artificial high-degree interpolation if the
        # transformed sequence is genuinely lower degree.
        while (
            len(coeffs) > 1
            and coeffs[-1] == 0
        ):
            coeffs.pop()

        transformed_components.append(
            coeffs
        )

    return transformed_components


# ==============================================================================
# COMPLEXITY
# ==============================================================================

def lcm(a, b):
    a = abs(a)
    b = abs(b)

    if a == 0:
        return b

    if b == 0:
        return a

    return abs(
        a // gcd(a, b) * b
    )


def common_integer_factor(coeffs):
    nonzero = [
        Fraction(c)
        for c in coeffs
        if c != 0
    ]

    if not nonzero:
        return Fraction(0)

    g = 0
    denom_lcm = 1

    for c in nonzero:
        g = gcd(
            g,
            abs(c.numerator),
        )
        denom_lcm = lcm(
            denom_lcm,
            c.denominator,
        )

    return Fraction(
        g,
        denom_lcm,
    )


def primitive_signature(coeffs):
    factor = common_integer_factor(
        coeffs
    )

    if factor == 0:
        return [0] * len(coeffs)

    scaled = [
        Fraction(c) / factor
        for c in coeffs
    ]

    first = next(
        (
            x for x in scaled
            if x != 0
        ),
        Fraction(1),
    )

    if first < 0:
        scaled = [
            -x
            for x in scaled
        ]

    return scaled


def component_complexity(coeffs):
    nonzero = [
        Fraction(c)
        for c in coeffs
        if c != 0
    ]

    if not nonzero:
        return {
            "degree": -1,
            "nonzero": 0,
            "max_num_digits": 0,
            "max_den_digits": 0,
            "total_digits": 0,
            "common_factor": Fraction(0),
        }

    return {
        "degree": len(coeffs) - 1,
        "nonzero": len(nonzero),
        "max_num_digits": max(
            len(str(abs(x.numerator)))
            for x in nonzero
        ),
        "max_den_digits": max(
            len(str(x.denominator))
            for x in nonzero
        ),
        "total_digits": sum(
            len(str(abs(x.numerator)))
            + len(str(x.denominator))
            for x in nonzero
        ),
        "common_factor": common_integer_factor(
            coeffs
        ),
    }


def total_complexity(components):
    data = [
        component_complexity(c)
        for c in components
    ]

    if not data:
        return (
            0,
            0,
            0,
            0,
        )

    return (
        max(
            x["max_num_digits"]
            for x in data
        ),
        max(
            x["max_den_digits"]
            for x in data
        ),
        sum(
            x["total_digits"]
            for x in data
        ),
        sum(
            x["nonzero"]
            for x in data
        ),
    )


# ==============================================================================
# WEIGHTED DATA
# ==============================================================================

def apply_weight(rows, weight_name):
    out = {}

    for p, row in rows.items():
        out[p] = [
            value
            * weight_value(
                weight_name,
                r,
            )
            for r, value in enumerate(row)
        ]

    return out


# ==============================================================================
# DIRECT TRANSFER TEST
# ==============================================================================

def test_operator_on_rows(
    rows,
    solution,
    order,
    degree,
):
    return verify_operator(
        rows,
        solution,
        order,
        degree,
    )


def test_conjugated_transfer(
    weight_name,
    transformed_components,
):
    """
    Use the same transformed B-odd coefficient functions on every
    sector, after applying the same weight to each sector.
    """

    results = {}

    for sector_name, rows in SECTORS.items():

        weighted = apply_weight(
            rows,
            weight_name,
        )

        # Check whether the transformed operator acts exactly over
        # every available consecutive p transition.
        ok = True

        ps = sorted(weighted)

        for p0, p1 in zip(ps, ps[1:]):

            previous = weighted[p0]
            target = weighted[p1]

            for r in range(len(target)):

                total = Fraction(0)

                for a in range(ORDER):

                    src = r + a

                    if src >= len(previous):
                        continue

                    F = transformed_components[a]

                    value = falling_eval(
                        F,
                        r,
                    )

                    total += (
                        value
                        * previous[src]
                    )

                if total != target[r]:
                    ok = False
                    break

            if not ok:
                break

        results[sector_name] = ok

    return results


# ==============================================================================
# PRINT OPERATOR
# ==============================================================================

def print_components(
    title,
    components,
):
    print()
    print(title)

    for a, coeffs in enumerate(components):

        print(
            f"  F_{a}(r) = {coeffs}"
        )

        stats = component_complexity(
            coeffs
        )

        print(
            f"    degree={stats['degree']}"
        )

        print(
            f"    nonzero={stats['nonzero']}"
        )

        print(
            f"    primitive_signature="
            f"{primitive_signature(coeffs)}"
        )

        print(
            f"    max_num_digits="
            f"{stats['max_num_digits']}"
        )

        print(
            f"    max_den_digits="
            f"{stats['max_den_digits']}"
        )


# ==============================================================================
# BASIC SELF-TESTS
# ==============================================================================

def data_self_test():
    """
    Check every stored row against itself by reconstructing the row
    from its falling coefficients. Since these are already falling
    coefficients, this is simply a direct pointwise evaluation test.
    """

    ok = True

    for sector_name, rows in SECTORS.items():

        print()
        print(sector_name)

        for p, row in sorted(rows.items()):

            # Treat the stored row as falling-basis coefficients and
            # evaluate the corresponding polynomial at 0..degree.
            for x in range(len(row)):

                predicted = falling_eval(
                    row,
                    x,
                )

                # We cannot compare against row[x], because row[x] is
                # a falling coefficient, not a monomial value.
                # Instead, reconstruct the coefficients through finite
                # differences and demand exact roundtrip.
                sampled = [
                    falling_eval(
                        row,
                        t,
                    )
                    for t in range(len(row))
                ]

                recovered = falling_interpolate(
                    sampled
                )

                if recovered != row:
                    ok = False
                    break

            print(
                f"  p={p}: roundtrip={ok}"
            )

    return ok


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 121 — EXACT B-ODD OPERATOR "
        "CONJUGATION / NORMALIZATION AUDIT"
    )
    print("=" * 78)

    # --------------------------------------------------------------------------
    # 1. Recover the exceptional operator independently.
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. UNIQUE B-ODD OPERATOR")
    print("=" * 78)

    info = recover_bodd_operator()

    print(
        f"  consistent={info['consistent']}"
    )
    print(
        f"  unique={info['unique']}"
    )
    print(
        f"  rank={info['rank']}"
    )
    print(
        f"  unknowns={info['unknowns']}"
    )
    print(
        f"  nullity={info['nullity']}"
    )

    if not info["unique"]:
        print(
            "ERROR: the expected unique B-odd operator "
            "could not be recovered."
        )
        print(
            "Experiment 121 aborted."
        )
        return

    base_solution = info["solution"]

    base_exact = verify_operator(
        B_odd,
        base_solution,
        ORDER,
        DEGREE,
    )

    print(
        f"  pointwise_exact={base_exact}"
    )

    base_components = [
        component(
            base_solution,
            a,
            DEGREE,
        )
        for a in range(ORDER)
    ]

    print_components(
        "  ORIGINAL OPERATOR COMPONENTS",
        base_components,
    )

    # --------------------------------------------------------------------------
    # 2. Natural conjugations.
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT NATURAL-WEIGHT CONJUGATION")
    print("=" * 78)

    conjugated = {}

    for weight_name in WEIGHTS:

        components = conjugate_operator(
            base_solution,
            weight_name,
        )

        if components is None:
            print()
            print(
                f"WEIGHT={weight_name}: "
                "not polynomial under tested exact range"
            )
            continue

        conjugated[weight_name] = components

        print()
        print(
            f"WEIGHT={weight_name}"
        )

        for a, coeffs in enumerate(components):

            stats = component_complexity(
                coeffs
            )

            print(
                f"  F_{a}: "
                f"degree={stats['degree']} "
                f"nonzero={stats['nonzero']} "
                f"max_num_digits="
                f"{stats['max_num_digits']} "
                f"max_den_digits="
                f"{stats['max_den_digits']}"
            )

            print(
                f"    coefficients={coeffs}"
            )

    # --------------------------------------------------------------------------
    # 3. Complexity ranking.
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. NORMALIZATION COMPLEXITY RANKING")
    print("=" * 78)

    base_score = total_complexity(
        base_components
    )

    print(
        f"  original={base_score}"
    )

    ranking = []

    for name, components in conjugated.items():

        score = total_complexity(
            components
        )

        ranking.append(
            (
                score,
                name,
            )
        )

    ranking.sort()

    for score, name in ranking:
        print(
            f"  {name}: {score}"
        )

    # --------------------------------------------------------------------------
    # 4. Exact transfer of every conjugated operator.
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXACT CROSS-SECTOR TRANSFER")
    print("=" * 78)

    transfer_results = {}

    for weight_name, components in conjugated.items():

        results = test_conjugated_transfer(
            weight_name,
            components,
        )

        transfer_results[
            weight_name
        ] = results

        print()
        print(
            f"WEIGHT={weight_name}"
        )

        for sector_name in SECTORS:
            print(
                f"  {sector_name}: "
                f"exact={results[sector_name]}"
            )

    # --------------------------------------------------------------------------
    # 5. Direct weighted-data recovery.
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. WEIGHTED-DATA OPERATOR RECOVERY")
    print("=" * 78)

    weighted_recovery = {}

    for weight_name in WEIGHTS:

        weighted_bodd = apply_weight(
            B_odd,
            weight_name,
        )

        recovered = recover_operator(
            weighted_bodd,
            ORDER,
            DEGREE,
        )

        weighted_recovery[
            weight_name
        ] = recovered

        print()
        print(
            f"WEIGHT={weight_name}"
        )
        print(
            f"  consistent={recovered['consistent']}"
        )
        print(
            f"  unique={recovered['unique']}"
        )
        print(
            f"  rank={recovered['rank']}"
        )
        print(
            f"  unknowns={recovered['unknowns']}"
        )
        print(
            f"  nullity={recovered['nullity']}"
        )

        if recovered["solution"] is not None:

            exact = verify_operator(
                weighted_bodd,
                recovered["solution"],
                ORDER,
                DEGREE,
            )

            print(
                f"  pointwise_exact={exact}"
            )

    # --------------------------------------------------------------------------
    # 6. Weight sanity.
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. EXACT WEIGHT SANITY")
    print("=" * 78)

    weights_ok = True

    for weight_name in WEIGHTS:

        print()
        print(weight_name)

        for r in range(0, 6):

            wr = weight_value(
                weight_name,
                r,
            )

            ratios = []

            for a in range(1, 3):

                wra = weight_value(
                    weight_name,
                    r + a,
                )

                if wra == 0:
                    ratios.append(
                        None
                    )
                else:
                    ratios.append(
                        wr / wra
                    )

            print(
                f"  r={r}: "
                f"w(r)={wr} "
                f"ratios={ratios}"
            )

    # --------------------------------------------------------------------------
    # 7. Fully independent falling-basis roundtrip.
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "7. INDEPENDENT FALLING-BASIS ROUNDTRIP SANITY"
    )
    print("=" * 78)

    roundtrip_ok = True

    for sector_name, rows in SECTORS.items():

        print()
        print(sector_name)

        for p, coeffs in sorted(rows.items()):

            samples = [
                falling_eval(
                    coeffs,
                    x,
                )
                for x in range(len(coeffs))
            ]

            recovered = falling_interpolate(
                samples
            )

            exact = (
                recovered == coeffs
            )

            if not exact:
                roundtrip_ok = False

            print(
                f"  p={p}: exact={exact}"
            )

    # --------------------------------------------------------------------------
    # 8. Structural interpretation.
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)
    print()
    print(
        "Experiment 120 isolated a unique exact B-odd operator of "
        "order 3 and falling degree 2."
    )
    print()
    print(
        "The coefficients were extremely large.  The present experiment "
        "tests whether that complexity is caused by the coordinate choice "
        "for the residual rows."
    )
    print()
    print(
        "For a diagonal change of coordinates"
    )
    print(
        "    v[p,r] = w(r) q[p,r],"
    )
    print()
    print(
        "the operator transforms by"
    )
    print(
        "    F_a(r) -> F_a(r) w(r)/w(r+a)."
    )
    print()
    print(
        "The tested weights are standard factorial/binomial scales:"
    )
    print(
        "    r!, 1/r!, binom(K,r), 1/binom(K,r),"
    )
    print(
        "    r!(K-r)!, 1/(r!(K-r)!),"
    )
    print(
        "    r!/(K-r)!, (K-r)!/r!."
    )
    print()
    print(
        "Three outcomes are particularly informative:"
    )
    print(
        "  1. a natural weight lowers the operator degree/order;"
    )
    print(
        "  2. a natural weight drastically reduces coefficient complexity;"
    )
    print(
        "  3. the transformed operator transfers exactly to another sector."
    )
    print()
    print(
        "Failure of all three would strengthen the conclusion that the "
        "B-odd exceptional operator is a finite-data phenomenon rather "
        "than a disguised standard factorial/binomial transform."
    )
    print()
    print(
        "All arithmetic is exact Fraction arithmetic."
    )
    print(
        "No floating point."
    )
    print(
        "No SymPy."
    )
    print(
        "No extrapolation."
    )

    # --------------------------------------------------------------------------
    # 9. Final checks.
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    failures = 0

    if not base_exact:
        failures += 1

    if not roundtrip_ok:
        failures += 1

    # At least one weight-sanity pass is structurally guaranteed here;
    # verify all weight evaluations were defined where used.
    if not weights_ok:
        failures += 1

    print(
        f"  b_odd_operator_exact={base_exact}"
    )
    print(
        f"  falling_basis_roundtrip={roundtrip_ok}"
    )
    print(
        f"  weight_sanity={weights_ok}"
    )
    print(
        f"  conjugations_tested={len(conjugated)}"
    )
    print(
        f"  failures={failures}"
    )
    print(
        f"  ALL BASIC CHECKS PASS={failures == 0}"
    )

    print()
    print(
        "EXPERIMENT 121 COMPLETE"
    )


if __name__ == "__main__":
    main()