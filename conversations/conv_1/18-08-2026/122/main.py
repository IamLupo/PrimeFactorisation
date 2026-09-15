#!/usr/bin/env python3

from fractions import Fraction
from math import gcd


# ==============================================================================
# EXPERIMENT 122R — EXACT B-ODD OPERATOR CONVENTION / FACTORIZATION AUDIT
#
# The supplied data are the corrected second-layer falling-basis rows
# from Experiment 115.
#
# We deliberately do NOT assume the shift convention.
#
# Tested forms include:
#
#   target[r] = sum_a F_a(r) source[r + a]
#   target[r] = sum_a F_a(r) source[r - a]
#
# and the reversed transition direction.
#
# F_a(r) is always reconstructed in the falling basis
#
#   F_a(r) = sum_d c[a,d] r_(d),
#
# with degree <= 2.
#
# A convention is accepted only when the exact system is consistent
# and reproduces EVERY available target entry.
# ==============================================================================


# ==============================================================================
# EXACT DATA
# ==============================================================================

ROWS = {
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


TRANSITIONS = [
    (1, 3),
    (3, 5),
    (5, 7),
]


# ==============================================================================
# FALLING FACTORIAL
# ==============================================================================

def falling(r, n):
    out = Fraction(1, 1)

    for i in range(n):
        out *= Fraction(r - i, 1)

    return out


# ==============================================================================
# POLYNOMIAL EVALUATION IN FALLING BASIS
# ==============================================================================

def eval_falling(coefficients, x):
    total = Fraction(0, 1)

    for d, coefficient in enumerate(coefficients):
        total += coefficient * falling(x, d)

    return total


# ==============================================================================
# FALLING-BASIS RECONSTRUCTION
# ==============================================================================

def values_to_falling(values):
    """
    Exact Newton/falling-factorial conversion:

        f(x) = sum_n c_n x_(n)

    where

        c_n = Delta^n f(0) / n!
    """

    work = list(values)
    coefficients = []
    factorial_value = 1

    while work:

        coefficients.append(
            work[0] / factorial_value
        )

        if len(work) == 1:
            break

        work = [
            work[i + 1] - work[i]
            for i in range(len(work) - 1)
        ]

        factorial_value *= len(coefficients)

    while (
        len(coefficients) > 1
        and coefficients[-1] == 0
    ):
        coefficients.pop()

    return coefficients


# ==============================================================================
# EXACT GAUSS-JORDAN
# ==============================================================================

def solve_exact(matrix, rhs):

    if not matrix:
        return {
            "consistent": True,
            "rank": 0,
            "unknowns": 0,
            "nullity": 0,
            "solution": [],
        }

    M = [
        list(map(Fraction, row)) + [Fraction(b)]
        for row, b in zip(matrix, rhs)
    ]

    nrows = len(M)
    ncols = len(M[0]) - 1

    pivot_columns = []
    row = 0

    for col in range(ncols):

        pivot = None

        for r in range(row, nrows):

            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        M[row], M[pivot] = M[pivot], M[row]

        pivot_value = M[row][col]

        M[row] = [
            x / pivot_value
            for x in M[row]
        ]

        for r in range(nrows):

            if r == row:
                continue

            factor = M[r][col]

            if factor == 0:
                continue

            M[r] = [
                M[r][j] - factor * M[row][j]
                for j in range(ncols + 1)
            ]

        pivot_columns.append(col)
        row += 1

        if row == nrows:
            break

    for r in range(nrows):

        left_zero = all(
            M[r][c] == 0
            for c in range(ncols)
        )

        if left_zero and M[r][-1] != 0:

            return {
                "consistent": False,
                "rank": len(pivot_columns),
                "unknowns": ncols,
                "nullity": None,
                "solution": None,
            }

    rank = len(pivot_columns)
    nullity = ncols - rank

    solution = None

    if nullity == 0:

        solution = [
            Fraction(0)
            for _ in range(ncols)
        ]

        for r, col in enumerate(pivot_columns):
            solution[col] = M[r][-1]

    return {
        "consistent": True,
        "rank": rank,
        "unknowns": ncols,
        "nullity": nullity,
        "solution": solution,
    }


# ==============================================================================
# OPERATOR APPLICATION
# ==============================================================================

def apply_operator(
    operator,
    source,
    r,
    shift_sign,
):
    """
    shift_sign = +1:

        source[r], source[r+1], source[r+2]

    shift_sign = -1:

        source[r], source[r-1], source[r-2]
    """

    total = Fraction(0, 1)

    for a in range(3):

        source_index = r + shift_sign * a

        if (
            source_index < 0
            or source_index >= len(source)
        ):
            continue

        total += (
            eval_falling(
                operator[a],
                r,
            )
            * source[source_index]
        )

    return total


# ==============================================================================
# OPERATOR SYSTEM
# ==============================================================================

def reconstruct_operator(
    transitions,
    shift_sign,
    reverse_direction,
    degree=2,
):
    """
    Unknowns:

        F_0(r): degree <= degree
        F_1(r): degree <= degree
        F_2(r): degree <= degree

    Total unknowns = 3*(degree+1).
    """

    matrix = []
    rhs = []

    for left_p, right_p in transitions:

        if reverse_direction:
            source_p = right_p
            target_p = left_p
        else:
            source_p = left_p
            target_p = right_p

        source = ROWS[source_p]
        target = ROWS[target_p]

        for r in range(len(target)):

            row = [
                Fraction(0)
                for _ in range(3 * (degree + 1))
            ]

            for a in range(3):

                source_index = (
                    r + shift_sign * a
                )

                if (
                    source_index < 0
                    or source_index >= len(source)
                ):
                    continue

                source_value = (
                    source[source_index]
                )

                for d in range(degree + 1):

                    column = (
                        a * (degree + 1) + d
                    )

                    row[column] += (
                        source_value
                        * falling(r, d)
                    )

            matrix.append(row)
            rhs.append(target[r])

    return solve_exact(matrix, rhs)


# ==============================================================================
# SPLIT OPERATOR
# ==============================================================================

def split_operator(solution, degree=2):

    width = degree + 1

    return [
        solution[
            0 * width:
            1 * width
        ],
        solution[
            1 * width:
            2 * width
        ],
        solution[
            2 * width:
            3 * width
        ],
    ]


# ==============================================================================
# POINTWISE VALIDATION
# ==============================================================================

def validate_operator(
    operator,
    shift_sign,
    reverse_direction,
):

    results = []

    for left_p, right_p in TRANSITIONS:

        if reverse_direction:
            source_p = right_p
            target_p = left_p
        else:
            source_p = left_p
            target_p = right_p

        source = ROWS[source_p]
        target = ROWS[target_p]

        exact = True

        for r in range(len(target)):

            predicted = apply_operator(
                operator,
                source,
                r,
                shift_sign,
            )

            if predicted != target[r]:
                exact = False
                break

        results.append(
            (
                source_p,
                target_p,
                exact,
            )
        )

    return results


# ==============================================================================
# POLYNOMIAL HELPERS
# ==============================================================================

def trim(poly):

    poly = list(poly)

    while (
        len(poly) > 1
        and poly[-1] == 0
    ):
        poly.pop()

    return poly


def poly_add(a, b):

    n = max(len(a), len(b))

    out = [
        Fraction(0)
        for _ in range(n)
    ]

    for i, x in enumerate(a):
        out[i] += x

    for i, x in enumerate(b):
        out[i] += x

    return trim(out)


def poly_multiply(a, b):

    if not a or not b:
        return [Fraction(0)]

    degree = (
        len(a) - 1
        + len(b) - 1
    )

    values = []

    for r in range(degree + 1):

        values.append(
            eval_falling(a, r)
            * eval_falling(b, r)
        )

    return values_to_falling(values)


def poly_shift(a, amount):

    degree = len(a) - 1

    values = [
        eval_falling(
            a,
            r + amount,
        )
        for r in range(degree + 1)
    ]

    return values_to_falling(values)


def poly_equal(a, b):

    n = max(len(a), len(b))

    for i in range(n):

        av = (
            a[i]
            if i < len(a)
            else Fraction(0)
        )

        bv = (
            b[i]
            if i < len(b)
            else Fraction(0)
        )

        if av != bv:
            return False

    return True


# ==============================================================================
# FIRST-ORDER FACTORIZATION
# ==============================================================================

def compose_first_order(
    A,
    B,
    C,
    D,
):
    """
    (A + B E)(C + D E)

      F_0 = A C

      F_1 = A D + B C(r+1)

      F_2 = B D(r+1)
    """

    C_shift = poly_shift(C, 1)
    D_shift = poly_shift(D, 1)

    F0 = poly_multiply(A, C)

    F1 = poly_add(
        poly_multiply(A, D),
        poly_multiply(B, C_shift),
    )

    F2 = poly_multiply(
        B,
        D_shift,
    )

    return [
        F0,
        F1,
        F2,
    ]


def factorization_exact(
    operator,
    A,
    B,
    C,
    D,
):

    candidate = compose_first_order(
        A,
        B,
        C,
        D,
    )

    return all(
        poly_equal(
            candidate[i],
            operator[i],
        )
        for i in range(3)
    )


# ==============================================================================
# SMALL FACTOR DICTIONARY
# ==============================================================================

def small_polynomials():

    return [
        [Fraction(-2)],
        [Fraction(-1)],
        [Fraction(1)],
        [Fraction(2)],

        [Fraction(0), Fraction(1)],
        [Fraction(1), Fraction(1)],
        [Fraction(-1), Fraction(1)],

        [
            Fraction(0),
            Fraction(0),
            Fraction(1),
        ],

        [
            Fraction(1),
            Fraction(0),
            Fraction(1),
        ],

        [
            Fraction(-1),
            Fraction(0),
            Fraction(1),
        ],
    ]


def search_factorizations(operator):

    candidates = []

    dictionary = small_polynomials()

    for A in dictionary:
        for B in dictionary:
            for C in dictionary:
                for D in dictionary:

                    if factorization_exact(
                        operator,
                        A,
                        B,
                        C,
                        D,
                    ):
                        candidates.append(
                            (A, B, C, D)
                        )

    return candidates


# ==============================================================================
# SIGNATURE
# ==============================================================================

def primitive_signature(poly):

    poly = trim(poly)

    denominator_lcm = 1

    for x in poly:

        d = x.denominator

        denominator_lcm = (
            denominator_lcm * d
            // gcd(denominator_lcm, d)
        )

    integers = [
        x.numerator
        * (
            denominator_lcm
            // x.denominator
        )
        for x in poly
    ]

    common = 0

    for x in integers:
        common = gcd(common, abs(x))

    if common == 0:
        return [0]

    integers = [
        x // common
        for x in integers
    ]

    for x in integers:

        if x != 0:

            if x < 0:
                integers = [
                    -y
                    for y in integers
                ]

            break

    return integers


# ==============================================================================
# PRINT OPERATOR
# ==============================================================================

def print_operator(operator):

    for i, F in enumerate(operator):

        print(
            f"  F_{i}: degree={len(F)-1}"
        )

        print(
            f"    coefficients={F}"
        )

        print(
            f"    signature={primitive_signature(F)}"
        )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 122R — EXACT B-ODD OPERATOR "
        "CONVENTION / FACTORIZATION AUDIT"
    )
    print("=" * 78)

    print()
    print("=" * 78)
    print("1. DATA VALIDATION")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        row = ROWS[p]

        print(
            f"  q_{p}: degree={len(row)-1} "
            f"entries={len(row)}"
        )

    # --------------------------------------------------------------------------
    # Convention search
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. OPERATOR-CONVENTION SEARCH")
    print("=" * 78)

    successful = []

    for reverse in [False, True]:

        for shift_sign in [1, -1]:

            result = reconstruct_operator(
                TRANSITIONS,
                shift_sign,
                reverse,
                degree=2,
            )

            label_direction = (
                "reverse"
                if reverse
                else "forward"
            )

            label_shift = (
                "r+a"
                if shift_sign == 1
                else "r-a"
            )

            print()
            print(
                f"  direction={label_direction} "
                f"shift={label_shift}"
            )

            print(
                f"    consistent={result['consistent']}"
            )

            print(
                f"    rank={result['rank']}"
            )

            print(
                f"    unknowns={result['unknowns']}"
            )

            print(
                f"    nullity={result['nullity']}"
            )

            if (
                result["consistent"]
                and result["solution"] is not None
                and result["nullity"] == 0
            ):

                operator = split_operator(
                    result["solution"]
                )

                checks = validate_operator(
                    operator,
                    shift_sign,
                    reverse,
                )

                exact = all(
                    ok
                    for _, _, ok in checks
                )

                print(
                    f"    pointwise_exact={exact}"
                )

                if exact:
                    successful.append(
                        (
                            reverse,
                            shift_sign,
                            operator,
                        )
                    )

    # --------------------------------------------------------------------------
    # Result
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. SELECTED EXACT CONVENTION")
    print("=" * 78)

    if not successful:

        print(
            "  NONE"
        )

        print()
        print(
            "  The supplied Experiment-115 rows do not admit a"
        )
        print(
            "  degree-2 three-term operator under any of the tested"
        )
        print(
            "  natural shift conventions."
        )

        print()
        print(
            "  Therefore the enormous coefficients printed in"
        )
        print(
            "  Experiments 118-121 cannot be safely reused as an"
        )
        print(
            "  operator on these rows without recovering the exact"
        )
        print(
            "  convention used by those scripts."
        )

        print()
        print("=" * 78)
        print("4. FINAL EXACTNESS")
        print("=" * 78)
        print(
            "  data_exact=True"
        )
        print(
            "  convention_search_completed=True"
        )
        print(
            "  exact_operator_found=False"
        )
        print(
            "  factorization_skipped=True"
        )
        print(
            "  ALL BASIC CHECKS PASS=False"
        )

        print()
        print(
            "EXPERIMENT 122R INCOMPLETE — OPERATOR CONVENTION UNRESOLVED"
        )

        return

    # There should normally be one.
    reverse, shift_sign, operator = successful[0]

    print(
        f"  direction={'reverse' if reverse else 'forward'}"
    )

    print(
        f"  shift={'r+a' if shift_sign == 1 else 'r-a'}"
    )

    print_operator(operator)

    # --------------------------------------------------------------------------
    # Factorization
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXACT FIRST-ORDER FACTORIZATION SEARCH")
    print("=" * 78)

    candidates = search_factorizations(
        operator
    )

    print(
        f"  exact_candidate_count={len(candidates)}"
    )

    for i, candidate in enumerate(
        candidates,
        start=1,
    ):

        A, B, C, D = candidate

        print(
            f"  candidate {i}:"
        )
        print(
            f"    A={A}"
        )
        print(
            f"    B={B}"
        )
        print(
            f"    C={C}"
        )
        print(
            f"    D={D}"
        )

    # --------------------------------------------------------------------------
    # Final
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  operator_reconstruction=True"
    )

    print(
        "  pointwise_exact=True"
    )

    print(
        "  factorization_search_completed=True"
    )

    print(
        f"  factorization_found={bool(candidates)}"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 122R COMPLETE"
    )


if __name__ == "__main__":
    main()

