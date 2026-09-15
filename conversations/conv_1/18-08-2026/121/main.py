#!/usr/bin/env python3

from fractions import Fraction
from math import factorial, gcd


# ==============================================================================
# EXPERIMENT 122 — EXACT B-ODD THREE-TERM SHIFT-OPERATOR
# RECONSTRUCTION / FACTORIZATION AUDIT
#
# Operator convention:
#
#     T = F_0(r) + F_1(r) E + F_2(r) E^2
#
# where
#
#     (E f)(r) = f(r+1).
#
# The operator acts on the exact residual falling-basis rows
#
#     q[p,r].
#
# For consecutive B-odd rows:
#
#     q_{p+2}(r)
#       = F_0(r) q_p(r)
#       + F_1(r) q_p(r+1)
#       + F_2(r) q_p(r+2).
#
# Everything is reconstructed from exact Fraction arithmetic.
# ==============================================================================


# ==============================================================================
# 1. EXACT B-ODD SECOND-LAYER ROWS FROM EXPERIMENT 115
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


# ==============================================================================
# 2. FALLING-FACTORIAL BASIS
# ==============================================================================

def falling(r, n):
    value = Fraction(1, 1)
    for i in range(n):
        value *= Fraction(r - i, 1)
    return value


def eval_falling(coeffs, r):
    total = Fraction(0, 1)

    for n, c in enumerate(coeffs):
        total += c * falling(r, n)

    return total


def values_to_falling(values):
    """
    Newton forward representation:

        f(r) = sum_n c_n r_(n)

    with

        c_n = Delta^n f(0) / n!
    """
    work = list(values)
    result = []

    while work:
        result.append(work[0])

        work = [
            work[i + 1] - work[i]
            for i in range(len(work) - 1)
        ]

        if not work:
            break

    # The values above are values in the ordinary integer basis.
    # Convert Newton forward differences to falling coefficients.
    #
    # Recompute independently to avoid accidental basis confusion.
    work = list(values)
    result = []

    n = 0
    while work:
        result.append(work[0] / factorial(n))

        work = [
            work[i + 1] - work[i]
            for i in range(len(work) - 1)
        ]

        n += 1

    while len(result) > 1 and result[-1] == 0:
        result.pop()

    return result


# ==============================================================================
# 3. EXACT LINEAR SYSTEM SOLVER
# ==============================================================================

def solve_exact(A, b):
    """
    Gauss-Jordan elimination over Fraction.

    Returns:
        consistent
        rank
        unknowns
        nullity
        solution, when unique
    """

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

    nrows = len(M)
    ncols = len(A[0])

    pivot_cols = []
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

        p = M[row][col]

        M[row] = [
            x / p
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

        pivot_cols.append(col)
        row += 1

        if row == nrows:
            break

    # Inconsistency check.
    for r in range(nrows):

        left_zero = all(
            M[r][c] == 0
            for c in range(ncols)
        )

        if left_zero and M[r][ncols] != 0:
            return {
                "consistent": False,
                "rank": len(pivot_cols),
                "unknowns": ncols,
                "nullity": None,
                "solution": None,
            }

    rank = len(pivot_cols)
    nullity = ncols - rank

    solution = None

    if nullity == 0:
        solution = [Fraction(0)] * ncols

        for i, col in enumerate(pivot_cols):
            solution[col] = M[i][ncols]

    return {
        "consistent": True,
        "rank": rank,
        "unknowns": ncols,
        "nullity": nullity,
        "solution": solution,
    }


# ==============================================================================
# 4. RECONSTRUCT THE THREE OPERATOR COMPONENTS
# ==============================================================================

def reconstruct_operator():

    degree = 2
    ncomp = 3

    # Unknowns:
    #
    # F_0(r) = a0 + a1 r_(1) + a2 r_(2)
    # F_1(r) = b0 + b1 r_(1) + b2 r_(2)
    # F_2(r) = c0 + c1 r_(1) + c2 r_(2)
    #
    # Nine unknowns total.

    A = []
    b = []

    transitions = [
        (1, 3),
        (3, 5),
        (5, 7),
    ]

    for source_p, target_p in transitions:

        source = ROWS[source_p]
        target = ROWS[target_p]

        # Target exists only through its displayed degree.
        for r in range(len(target)):

            row = [Fraction(0)] * 9

            for a in range(ncomp):

                q_index = r + a

                if q_index >= len(source):
                    continue

                q_value = source[q_index]

                for d in range(degree + 1):

                    row[
                        a * (degree + 1) + d
                    ] += (
                        q_value
                        * falling(r, d)
                    )

            A.append(row)
            b.append(target[r])

    return solve_exact(A, b)


def split_operator(solution):

    F0 = solution[0:3]
    F1 = solution[3:6]
    F2 = solution[6:9]

    return [F0, F1, F2]


# ==============================================================================
# 5. EXACT OPERATOR EVALUATION
# ==============================================================================

def operator_value(operator, source_row, r):

    total = Fraction(0, 1)

    for shift, F in enumerate(operator):

        idx = r + shift

        if idx >= len(source_row):
            continue

        total += (
            eval_falling(F, r)
            * source_row[idx]
        )

    return total


def pointwise_check(operator):

    checks = []

    for source_p, target_p in [
        (1, 3),
        (3, 5),
        (5, 7),
    ]:

        source = ROWS[source_p]
        target = ROWS[target_p]

        ok = True

        for r in range(len(target)):

            predicted = operator_value(
                operator,
                source,
                r,
            )

            if predicted != target[r]:
                ok = False
                break

        checks.append(
            (source_p, target_p, ok)
        )

    return checks


# ==============================================================================
# 6. EXACT POLYNOMIAL OPERATIONS
# ==============================================================================

def trim(poly):

    poly = list(poly)

    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()

    if not poly:
        return [Fraction(0)]

    return poly


def add_poly(a, b):

    n = max(len(a), len(b))

    out = [Fraction(0)] * n

    for i, x in enumerate(a):
        out[i] += x

    for i, x in enumerate(b):
        out[i] += x

    return trim(out)


def multiply_poly(a, b):

    da = len(a) - 1
    db = len(b) - 1

    values = []

    for r in range(da + db + 1):

        ar = eval_falling(a, r)
        br = eval_falling(b, r)

        values.append(ar * br)

    return values_to_falling(values)


def shift_poly(a, shift):

    degree = len(a) - 1

    values = [
        eval_falling(a, r + shift)
        for r in range(degree + 1)
    ]

    return values_to_falling(values)


def poly_equal(a, b):

    n = max(len(a), len(b))

    for i in range(n):

        ai = a[i] if i < len(a) else Fraction(0)
        bi = b[i] if i < len(b) else Fraction(0)

        if ai != bi:
            return False

    return True


# ==============================================================================
# 7. FIRST-ORDER SHIFT COMPOSITION
# ==============================================================================

def compose_first_order(A, B, C, D):
    """
    (A + B E)(C + D E)

      = AC
      + [AD + B*C(r+1)] E
      + [B*D(r+1)] E^2.
    """

    C1 = shift_poly(C, 1)
    D1 = shift_poly(D, 1)

    G0 = multiply_poly(A, C)

    G1 = add_poly(
        multiply_poly(A, D),
        multiply_poly(B, C1),
    )

    G2 = multiply_poly(
        B,
        D1,
    )

    return [G0, G1, G2]


# ==============================================================================
# 8. EXACT FIRST-ORDER FACTORIZATION TEST
# ==============================================================================

def complete_factorization_check(
    A, B, C, D,
    target,
):

    candidate = compose_first_order(
        A, B, C, D
    )

    return all(
        poly_equal(candidate[i], target[i])
        for i in range(3)
    )


# ==============================================================================
# 9. SMALL FACTOR DICTIONARY
# ==============================================================================

def factor_dictionary():

    factors = []

    # Constants
    for c in [
        Fraction(-2),
        Fraction(-1),
        Fraction(1),
        Fraction(2),
    ]:
        factors.append([c])

    # degree 1
    factors.extend([
        [Fraction(0), Fraction(1)],      # r
        [Fraction(1), Fraction(1)],      # r+1
        [Fraction(-1), Fraction(1)],     # r-1
    ])

    # degree 2
    factors.extend([
        [Fraction(0), Fraction(0), Fraction(1)],   # r_(2)
        [Fraction(1), Fraction(0), Fraction(1)],
        [Fraction(-1), Fraction(0), Fraction(1)],
    ])

    return factors


# ==============================================================================
# 10. SEARCH FACTORIZATION BY EXPLICIT COMPOSITION
# ==============================================================================

def search_small_factorizations(target):

    candidates = []

    dictionary = factor_dictionary()

    # Test all combinations with A,C chosen from the small dictionary,
    # and B,D chosen likewise.
    #
    # This is deliberately finite.

    for A in dictionary:
        for B in dictionary:
            for C in dictionary:
                for D in dictionary:

                    # Avoid completely zero first-order operators.
                    if all(x == 0 for x in A) and all(
                        x == 0 for x in B
                    ):
                        continue

                    if all(x == 0 for x in C) and all(
                        x == 0 for x in D
                    ):
                        continue

                    if complete_factorization_check(
                        A, B, C, D, target
                    ):
                        candidates.append(
                            (A, B, C, D)
                        )

    return candidates


# ==============================================================================
# 11. EXACT PRIMITIVE SIGNATURE
# ==============================================================================

def primitive_signature(poly):

    poly = trim(poly)

    common_den = 1

    for x in poly:
        common_den = (
            common_den * x.denominator
            // gcd(common_den, x.denominator)
        )

    integers = [
        x.numerator
        * (common_den // x.denominator)
        for x in poly
    ]

    g = 0

    for x in integers:
        g = gcd(g, abs(x))

    if g == 0:
        return [0]

    integers = [
        x // g
        for x in integers
    ]

    # Canonical sign.
    for x in integers:
        if x != 0:
            if x < 0:
                integers = [-y for y in integers]
            break

    return integers


# ==============================================================================
# 12. MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 122 — EXACT B-ODD THREE-TERM "
        "SHIFT-OPERATOR FACTORIZATION AUDIT"
    )
    print("=" * 78)

    # --------------------------------------------------------------------------
    # Reconstruction
    # --------------------------------------------------------------------------

    result = reconstruct_operator()

    print()
    print("=" * 78)
    print("1. EXACT OPERATOR RECONSTRUCTION")
    print("=" * 78)

    print(
        f"  consistent={result['consistent']}"
    )
    print(
        f"  rank={result['rank']}"
    )
    print(
        f"  unknowns={result['unknowns']}"
    )
    print(
        f"  nullity={result['nullity']}"
    )

    if not result["consistent"] or result["solution"] is None:
        print()
        print("  OPERATOR RECONSTRUCTION FAILED")
        return

    operator = split_operator(
        result["solution"]
    )

    for i, F in enumerate(operator):

        print(
            f"  F_{i}: degree={len(F)-1}"
        )

        print(
            f"    coefficients={F}"
        )

        print(
            f"    primitive_signature={primitive_signature(F)}"
        )

    # --------------------------------------------------------------------------
    # Pointwise verification
    # --------------------------------------------------------------------------

    checks = pointwise_check(
        operator
    )

    print()
    print("=" * 78)
    print("2. POINTWISE OPERATOR RECONSTRUCTION")
    print("=" * 78)

    for source_p, target_p, ok in checks:

        print(
            f"  {source_p}->{target_p}: exact={ok}"
        )

    operator_exact = all(
        ok for _, _, ok in checks
    )

    # --------------------------------------------------------------------------
    # Factorization
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. COMPLETE FIRST-ORDER FACTORIZATION SEARCH")
    print("=" * 78)

    print()
    print(
        "Testing exact identities"
    )
    print()
    print(
        "  F_0 = A*C"
    )
    print(
        "  F_1 = A*D + B*C(r+1)"
    )
    print(
        "  F_2 = B*D(r+1)"
    )

    candidates = search_small_factorizations(
        operator
    )

    print()
    print(
        f"  exact_small_factorization_count={len(candidates)}"
    )

    for i, (A, B, C, D) in enumerate(candidates):

        print(
            f"  candidate {i+1}:"
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

    if not candidates:
        print(
            "  NONE"
        )

    # --------------------------------------------------------------------------
    # Trivial factorization audit
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. TRIVIAL IDENTITY FACTORIZATION AUDIT")
    print("=" * 78)

    identity = []

    zero = [Fraction(0)]
    one = [Fraction(1)]

    # T = T * I
    if complete_factorization_check(
        operator[0],
        operator[1],
        one,
        zero,
        operator,
    ):
        identity.append(
            "T = T * I"
        )

    # T = I * T
    if complete_factorization_check(
        one,
        zero,
        operator[0],
        operator[1],
        operator,
    ):
        identity.append(
            "T = I * T"
        )

    print(
        f"  identities={identity}"
    )

    # --------------------------------------------------------------------------
    # Structural interpretation
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print()
    print(
        "The previous run failed because the operator coefficients were"
    )
    print(
        "not reconstructed from the exact B-odd data before the"
    )
    print(
        "factorization test."
    )

    print()
    print(
        "This experiment reconstructs the operator independently from"
    )
    print(
        "the four exact B-odd residual rows:"
    )

    print()
    print(
        "    q_3 = T[q_1]"
    )
    print(
        "    q_5 = T[q_3]"
    )
    print(
        "    q_7 = T[q_5]"
    )

    print()
    print(
        "Only after that reconstruction does it test"
    )
    print()
    print(
        "    T = (A + B E)(C + D E)"
    )
    print()
    print(
        "using the complete three coefficient identities."
    )

    print()
    print(
        "A candidate is therefore reported only when the entire"
    )
    print(
        "operator is reproduced exactly."
    )

    print()
    print(
        "The finite factor dictionary is intentionally small."
    )
    print(
        "Failure of this search is not a proof of irreducibility"
    )
    print(
        "over the full rational Ore-polynomial ring."
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
    # Final exactness
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    failures = 0

    if not operator_exact:
        failures += 1

    print(
        f"  operator_reconstruction={operator_exact}"
    )

    print(
        "  factorization_audit_completed=True"
    )

    print(
        f"  exact_small_factorization_found={bool(candidates)}"
    )

    print(
        f"  failures={failures}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={failures == 0}"
    )

    print()
    print(
        "EXPERIMENT 122 COMPLETE"
    )


if __name__ == "__main__":
    main()

