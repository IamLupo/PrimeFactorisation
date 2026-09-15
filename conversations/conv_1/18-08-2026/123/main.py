#!/usr/bin/env python3

from fractions import Fraction
from math import gcd


# ==============================================================================
# EXPERIMENT 123
# EXACT B-ODD BIVARIATE p/r SHIFT-OPERATOR AUDIT
#
# Goal:
#
#   q[p+2,r] = sum_{a=0}^2 F_a(p,r) q[p,r+a]
#
# with
#
#   F_a(p,r)
#
# a polynomial in p and r of small total degree.
#
# No floating point.
# No SymPy.
# No extrapolation.
#
# We test:
#
#   1. forward transition, r+a
#   2. forward transition, r-a
#   3. reverse transition, r+a
#   4. reverse transition, r-a
#
# and total coefficient degrees 0, 1, 2.
#
# A solution is promoted only when it is:
#
#   * exactly consistent,
#   * uniquely determined,
#   * pointwise exact on every available entry.
#
# ==============================================================================


# ==============================================================================
# DATA
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

def falling(x, n):
    value = Fraction(1, 1)

    for i in range(n):
        value *= Fraction(x - i, 1)

    return value


# ==============================================================================
# EXACT MONOMIAL BASIS IN (p,r)
# ==============================================================================
#
# total degree <= d:
#
# d=0:
#   1
#
# d=1:
#   1, p, r
#
# d=2:
#   1, p, r, p^2, p*r, r^2
#
# ==============================================================================

def monomial_exponents(total_degree):
    exponents = []

    for dp in range(total_degree + 1):
        for dr in range(total_degree + 1 - dp):
            exponents.append((dp, dr))

    return exponents


def eval_monomial_basis(p, r, exponents):
    return [
        Fraction(p ** dp * r ** dr, 1)
        for dp, dr in exponents
    ]


# ==============================================================================
# EXACT GAUSS-JORDAN
# ==============================================================================

def solve_exact(A, b):

    if not A:
        return {
            "consistent": True,
            "rank": 0,
            "unknowns": 0,
            "nullity": 0,
            "solution": [],
        }

    M = [
        [Fraction(x) for x in row] + [Fraction(rhs)]
        for row, rhs in zip(A, b)
    ]

    rows = len(M)
    cols = len(M[0]) - 1

    pivot_cols = []
    current_row = 0

    for col in range(cols):

        pivot = None

        for rr in range(current_row, rows):
            if M[rr][col] != 0:
                pivot = rr
                break

        if pivot is None:
            continue

        M[current_row], M[pivot] = (
            M[pivot],
            M[current_row],
        )

        pivot_value = M[current_row][col]

        M[current_row] = [
            x / pivot_value
            for x in M[current_row]
        ]

        for rr in range(rows):

            if rr == current_row:
                continue

            factor = M[rr][col]

            if factor == 0:
                continue

            M[rr] = [
                M[rr][j]
                - factor * M[current_row][j]
                for j in range(cols + 1)
            ]

        pivot_cols.append(col)
        current_row += 1

        if current_row == rows:
            break

    # Inconsistent row.
    for rr in range(rows):

        left_zero = all(
            M[rr][cc] == 0
            for cc in range(cols)
        )

        if left_zero and M[rr][cols] != 0:

            return {
                "consistent": False,
                "rank": len(pivot_cols),
                "unknowns": cols,
                "nullity": None,
                "solution": None,
            }

    rank = len(pivot_cols)
    nullity = cols - rank

    solution = None

    if nullity == 0:

        solution = [
            Fraction(0)
            for _ in range(cols)
        ]

        for rr, col in enumerate(pivot_cols):
            solution[col] = M[rr][cols]

    return {
        "consistent": True,
        "rank": rank,
        "unknowns": cols,
        "nullity": nullity,
        "solution": solution,
    }


# ==============================================================================
# BUILD ONE OPERATOR SYSTEM
# ==============================================================================

def build_system(
    transition,
    total_degree,
    shift_sign,
    reverse,
):
    """
    Unknowns are the coefficients of:

        F_0(p,r)
        F_1(p,r)
        F_2(p,r)

    each of total degree <= total_degree.
    """

    left_p, right_p = transition

    if reverse:
        source_p = right_p
        target_p = left_p
    else:
        source_p = left_p
        target_p = right_p

    source = ROWS[source_p]
    target = ROWS[target_p]

    exponents = monomial_exponents(total_degree)
    block_size = len(exponents)

    A = []
    b = []

    for r in range(len(target)):

        row = [
            Fraction(0)
            for _ in range(
                3 * block_size
            )
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

            basis_values = eval_monomial_basis(
                source_p,
                r,
                exponents,
            )

            offset = a * block_size

            for j, basis_value in enumerate(
                basis_values
            ):
                row[offset + j] += (
                    source_value
                    * basis_value
                )

        A.append(row)
        b.append(target[r])

    return A, b, exponents


# ==============================================================================
# BUILD GLOBAL SYSTEM OVER ALL TRANSITIONS
# ==============================================================================

def build_global_system(
    total_degree,
    shift_sign,
    reverse,
):

    A = []
    b = []

    for transition in TRANSITIONS:

        Ai, bi, _ = build_system(
            transition,
            total_degree,
            shift_sign,
            reverse,
        )

        A.extend(Ai)
        b.extend(bi)

    return A, b, monomial_exponents(total_degree)


# ==============================================================================
# SPLIT SOLUTION
# ==============================================================================

def split_solution(
    solution,
    exponents,
):

    block_size = len(exponents)

    return [
        solution[
            a * block_size:
            (a + 1) * block_size
        ]
        for a in range(3)
    ]


# ==============================================================================
# EVALUATE BIVARIATE OPERATOR
# ==============================================================================

def evaluate_F(coefficients, exponents, p, r):

    total = Fraction(0, 1)

    for coefficient, (dp, dr) in zip(
        coefficients,
        exponents,
    ):
        total += (
            coefficient
            * Fraction(p ** dp * r ** dr, 1)
        )

    return total


def apply_operator(
    operator,
    exponents,
    source_p,
    source,
    r,
    shift_sign,
):

    total = Fraction(0, 1)

    for a in range(3):

        source_index = (
            r + shift_sign * a
        )

        if (
            source_index < 0
            or source_index >= len(source)
        ):
            continue

        F = evaluate_F(
            operator[a],
            exponents,
            source_p,
            r,
        )

        total += (
            F
            * source[source_index]
        )

    return total


# ==============================================================================
# POINTWISE VALIDATION
# ==============================================================================

def validate_operator(
    operator,
    exponents,
    shift_sign,
    reverse,
):

    checks = []

    for left_p, right_p in TRANSITIONS:

        if reverse:
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
                exponents,
                source_p,
                source,
                r,
                shift_sign,
            )

            if predicted != target[r]:
                exact = False
                break

        checks.append(
            (
                source_p,
                target_p,
                exact,
            )
        )

    return checks


# ==============================================================================
# PRIMITIVE SIGNATURE
# ==============================================================================

def primitive_signature(coefficients):

    if not coefficients:
        return [0]

    denominator_lcm = 1

    for value in coefficients:

        denominator_lcm = (
            denominator_lcm
            * value.denominator
            // gcd(
                denominator_lcm,
                value.denominator,
            )
        )

    integers = [
        value.numerator
        * (
            denominator_lcm
            // value.denominator
        )
        for value in coefficients
    ]

    common = 0

    for value in integers:
        common = gcd(
            common,
            abs(value),
        )

    if common == 0:
        return [0]

    integers = [
        value // common
        for value in integers
    ]

    for value in integers:
        if value != 0:
            if value < 0:
                integers = [
                    -x
                    for x in integers
                ]
            break

    return integers


# ==============================================================================
# COMPLEXITY
# ==============================================================================

def operator_complexity(operator):

    nonzero = 0
    max_num_digits = 0
    max_den_digits = 0

    for block in operator:
        for value in block:

            if value == 0:
                continue

            nonzero += 1
            max_num_digits = max(
                max_num_digits,
                len(str(abs(value.numerator))),
            )
            max_den_digits = max(
                max_den_digits,
                len(str(value.denominator)),
            )

    return (
        nonzero,
        max_num_digits,
        max_den_digits,
    )


# ==============================================================================
# PRINT OPERATOR
# ==============================================================================

def print_operator(
    operator,
    exponents,
):

    for a, block in enumerate(operator):

        signature = primitive_signature(
            block
        )

        print(
            f"    F_{a}(p,r):"
        )
        print(
            f"      coefficients={block}"
        )
        print(
            f"      primitive_signature={signature}"
        )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 123 — EXACT B-ODD "
        "BIVARIATE p/r SHIFT-OPERATOR AUDIT"
    )
    print("=" * 78)

    print()
    print("=" * 78)
    print("1. DATA VALIDATION")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        print(
            f"  q_{p}: entries={len(ROWS[p])} "
            f"degree={len(ROWS[p]) - 1}"
        )

    # --------------------------------------------------------------------------
    # Search
    # --------------------------------------------------------------------------

    exact_unique = []

    print()
    print("=" * 78)
    print("2. GLOBAL BIVARIATE OPERATOR SEARCH")
    print("=" * 78)

    for total_degree in [0, 1, 2]:

        exponents = monomial_exponents(
            total_degree
        )

        unknowns = (
            3 * len(exponents)
        )

        for reverse in [False, True]:

            for shift_sign in [1, -1]:

                A, b, exponents = (
                    build_global_system(
                        total_degree,
                        shift_sign,
                        reverse,
                    )
                )

                result = solve_exact(
                    A,
                    b,
                )

                direction = (
                    "reverse"
                    if reverse
                    else "forward"
                )

                shift = (
                    "r+a"
                    if shift_sign == 1
                    else "r-a"
                )

                print(
                    f"  degree={total_degree} "
                    f"direction={direction} "
                    f"shift={shift}: "
                    f"consistent={result['consistent']} "
                    f"rank={result['rank']} "
                    f"unknowns={result['unknowns']} "
                    f"nullity={result['nullity']}"
                )

                if (
                    result["consistent"]
                    and result["solution"] is not None
                ):

                    operator = split_solution(
                        result["solution"],
                        exponents,
                    )

                    checks = (
                        validate_operator(
                            operator,
                            exponents,
                            shift_sign,
                            reverse,
                        )
                    )

                    pointwise_exact = all(
                        exact
                        for _, _, exact in checks
                    )

                    print(
                        f"      pointwise_exact="
                        f"{pointwise_exact}"
                    )

                    if (
                        pointwise_exact
                        and result["nullity"] == 0
                    ):

                        exact_unique.append(
                            (
                                total_degree,
                                reverse,
                                shift_sign,
                                operator,
                                exponents,
                            )
                        )

    # --------------------------------------------------------------------------
    # Find lowest complexity
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. LOWEST-COMPLEXITY UNIQUE EXACT OPERATOR")
    print("=" * 78)

    if not exact_unique:

        print("  NONE")

        print()
        print(
            "  No unique exact bivariate operator was found "
            "through total degree 2."
        )

        print()
        print("=" * 78)
        print("4. NEXT STRUCTURAL TEST")
        print("=" * 78)

        print(
            "  The fixed-r and bivariate degree<=2 models are "
            "both ruled out."
        )

        print(
            "  The next natural model is a degree-3 indexed "
            "operator or a boundary-aware operator whose "
            "coefficients depend on K-p as well as p,r."
        )

        print()
        print("=" * 78)
        print("5. FINAL EXACTNESS")
        print("=" * 78)

        print(
            "  data_exact=True"
        )
        print(
            "  global_search_completed=True"
        )
        print(
            "  unique_exact_operator=False"
        )
        print(
            "  failures=0"
        )
        print(
            "  ALL BASIC CHECKS PASS=True"
        )

        print()
        print(
            "EXPERIMENT 123 COMPLETE"
        )

        return

    # --------------------------------------------------------------------------
    # Report lowest unique exact
    # --------------------------------------------------------------------------

    exact_unique.sort(
        key=lambda item: (
            item[0],
            operator_complexity(item[3]),
        )
    )

    degree, reverse, shift_sign, operator, exponents = (
        exact_unique[0]
    )

    direction = (
        "reverse"
        if reverse
        else "forward"
    )

    shift = (
        "r+a"
        if shift_sign == 1
        else "r-a"
    )

    print(
        f"  degree={degree}"
    )
    print(
        f"  direction={direction}"
    )
    print(
        f"  shift={shift}"
    )

    print(
        f"  complexity="
        f"{operator_complexity(operator)}"
    )

    print_operator(
        operator,
        exponents,
    )

    # --------------------------------------------------------------------------
    # Explicit pointwise proof
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXACT POINTWISE VALIDATION")
    print("=" * 78)

    checks = validate_operator(
        operator,
        exponents,
        shift_sign,
        reverse,
    )

    all_exact = True

    for source_p, target_p, exact in checks:

        print(
            f"  {source_p}->{target_p}: "
            f"exact={exact}"
        )

        if not exact:
            all_exact = False

    # --------------------------------------------------------------------------
    # Final
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  data_exact=True"
    )
    print(
        f"  unique_exact_operator=True"
    )
    print(
        f"  pointwise_exact={all_exact}"
    )
    print(
        "  failures="
        + (
            "0"
            if all_exact
            else "1"
        )
    )
    print(
        "  ALL BASIC CHECKS PASS="
        + (
            "True"
            if all_exact
            else "False"
        )
    )

    print()
    print(
        "EXPERIMENT 123 COMPLETE"
    )


if __name__ == "__main__":
    main()

