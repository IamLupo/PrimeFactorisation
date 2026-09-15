#!/usr/bin/env python3

from fractions import Fraction
from math import gcd, log10
from functools import reduce


# =============================================================================
# EXPERIMENT 131 — EXACT B-ODD OPERATOR DETERMINANT / CRAMER-SCALE AUDIT
# =============================================================================
#
# PURPOSE
# -------
# Fact-check the hypothesis that the enormous rational coefficients in the
# exact boundary-defect operator are primarily caused by solving a full-rank
# 14 x 14 rational linear system.
#
# We reconstruct exactly:
#
#     q_{p+2}(r)
#       =
#       F_0(p,d) q_p(r)
#       +
#       F_1(p,d) q_p(r+1)
#       +
#       H_0(p) [d=0],
#
# where
#
#     d = D(p) - r.
#
# The ansatz is
#
#     F_a(p,d) = sum_{i+j<=2} c_{a,i,j} p^i d^j
#
# giving 6 coefficients for F_0,
#       6 coefficients for F_1,
#       2 coefficients for H_0,
#
# for a total of 14 unknowns.
#
# The exact finite data give:
#
#     p=1 -> p=3 : 6 equations
#     p=3 -> p=5 : 5 equations
#     p=5 -> p=7 : 3 equations
#
# total = 14 equations.
#
# THIS EXPERIMENT DOES NOT LOAD ANY PROJECT FILE.
#
# Everything is embedded directly below.
#
# No SymPy.
# No floating point in the mathematics.
# No extrapolation.
# =============================================================================


# =============================================================================
# 1. EXACT DATA
# =============================================================================

Q = {
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


D = {
    1: 5,
    3: 4,
    5: 2,
}


# =============================================================================
# 2. BASIS
# =============================================================================

# Total degree <= 2 in (p,d):
#
# 1, d, d^2, p, p*d, p^2
#
# Keep this order fixed throughout the determinant audit.

BASIS = [
    (0, 0),
    (0, 1),
    (0, 2),
    (1, 0),
    (1, 1),
    (2, 0),
]


# Unknown order:
#
# F0: six coefficients
# F1: six coefficients
# H0: two coefficients

UNKNOWN_NAMES = []

for i, j in BASIS:
    UNKNOWN_NAMES.append("F0_p%d_d%d" % (i, j))

for i, j in BASIS:
    UNKNOWN_NAMES.append("F1_p%d_d%d" % (i, j))

UNKNOWN_NAMES.append("H0_1")
UNKNOWN_NAMES.append("H0_p")


# =============================================================================
# 3. BASIC EXACT UTILITIES
# =============================================================================

def falling_values():
    return all(
        isinstance(value, Fraction)
        for row in Q.values()
        for value in row
    )


def q_value(p, r):
    if r < 0:
        return Fraction(0)
    row = Q[p]
    if r >= len(row):
        return Fraction(0)
    return row[r]


def monomial_basis_value(p, d):
    return [
        Fraction(1),
        Fraction(d),
        Fraction(d * d),
        Fraction(p),
        Fraction(p * d),
        Fraction(p * p),
    ]


def is_zero_row(row):
    return all(x == 0 for x in row)


# =============================================================================
# 4. BUILD THE EXACT 14 x 14 SYSTEM
# =============================================================================

def build_system():
    rows = []
    rhs = []
    meta = []

    transitions = [
        (1, 3),
        (3, 5),
        (5, 7),
    ]

    for p, p_next in transitions:
        terminal = D[p]

        # Include every available r in q_p.
        #
        # This automatically gives:
        #
        # 1 -> 3 : r = 0,...,5  -> 6 equations
        # 3 -> 5 : r = 0,...,4  -> 5 equations
        # 5 -> 7 : r = 0,...,2  -> 3 equations
        #
        # total = 14.

        for r in range(len(Q[p])):
            d = terminal - r
            basis = monomial_basis_value(p, d)

            row = []

            # F0 block
            for value in basis:
                row.append(value * q_value(p, r))

            # F1 block
            for value in basis:
                row.append(value * q_value(p, r + 1))

            # H0 block
            boundary = 1 if d == 0 else 0
            row.append(Fraction(boundary))
            row.append(Fraction(p * boundary))

            target = q_value(p_next, r)

            rows.append(row)
            rhs.append(target)
            meta.append((p, p_next, r, d))

    return rows, rhs, meta


# =============================================================================
# 5. ROW-CLEARING
# =============================================================================

def row_clear_fraction_matrix(A, b):
    """
    Convert a rational matrix [A|b] into an integer matrix by clearing each
    row denominator separately.

    Returns:
        integer_matrix,
        row_scale_product

    so that

        det(A) = det(integer_matrix) / row_scale_product

    exactly.
    """
    integer_rows = []
    scales = []

    for row, target in zip(A, b):
        denoms = [x.denominator for x in row]
        denoms.append(target.denominator)

        scale = 1
        for d in denoms:
            scale = lcm(scale, d)

        integer_row = [int(x * scale) for x in row]
        integer_rows.append(integer_row)
        scales.append(scale)

    product_scale = 1
    for s in scales:
        product_scale *= s

    return integer_rows, rhs_clear(b, scales), product_scale


def rhs_clear(b, scales):
    return [
        int(value * scale)
        for value, scale in zip(b, scales)
    ]


def lcm(a, b):
    if a == 0 or b == 0:
        return 0
    return abs(a // gcd(a, b) * b)


# =============================================================================
# 6. EXACT INTEGER DETERMINANT — BAREISS
# =============================================================================

def bareiss_det(A):
    """
    Fraction-free exact determinant for an integer matrix.
    """
    n = len(A)

    if n == 0:
        return 1

    if n == 1:
        return A[0][0]

    M = [row[:] for row in A]
    sign = 1
    prev = 1

    for k in range(n - 1):

        pivot = k
        while pivot < n and M[pivot][k] == 0:
            pivot += 1

        if pivot == n:
            return 0

        if pivot != k:
            M[k], M[pivot] = M[pivot], M[k]
            sign *= -1

        pivot_value = M[k][k]

        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = (
                    M[i][j] * pivot_value
                    - M[i][k] * M[k][j]
                )

                if k == 0:
                    M[i][j] = numerator
                else:
                    if numerator % prev != 0:
                        raise ArithmeticError(
                            "Bareiss exact division failed."
                        )
                    M[i][j] = numerator // prev

        for i in range(k + 1, n):
            M[i][k] = 0

        prev = pivot_value

    return sign * M[n - 1][n - 1]


# =============================================================================
# 7. EXACT RATIONAL GAUSSIAN SOLVER
# =============================================================================

def solve_exact(A, b):
    A = [[Fraction(x) for x in row] for row in A]
    b = [Fraction(x) for x in b]

    n = len(A)
    m = len(A[0])

    M = [
        A[i] + [b[i]]
        for i in range(n)
    ]

    row = 0
    pivots = []

    for col in range(m):
        pivot = None

        for r in range(row, n):
            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        M[row], M[pivot] = M[pivot], M[row]

        pv = M[row][col]

        for j in range(col, m + 1):
            M[row][j] /= pv

        for r in range(n):
            if r == row:
                continue

            factor = M[r][col]

            if factor == 0:
                continue

            for j in range(col, m + 1):
                M[r][j] -= factor * M[row][j]

        pivots.append(col)
        row += 1

    # Inconsistency
    for r in range(n):
        if all(M[r][c] == 0 for c in range(m)):
            if M[r][m] != 0:
                return None

    if len(pivots) != m:
        return None

    solution = [Fraction(0) for _ in range(m)]

    for r, c in enumerate(pivots):
        solution[c] = M[r][m]

    return solution


# =============================================================================
# 8. MATRIX RANK
# =============================================================================

def rank_exact(A):
    if not A:
        return 0

    M = [[Fraction(x) for x in row] for row in A]

    rows = len(M)
    cols = len(M[0])

    rank = 0

    for col in range(cols):
        pivot = None

        for r in range(rank, rows):
            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        M[rank], M[pivot] = M[pivot], M[rank]

        pv = M[rank][col]

        for j in range(col, cols):
            M[rank][j] /= pv

        for r in range(rows):
            if r == rank:
                continue

            factor = M[r][col]

            if factor == 0:
                continue

            for j in range(col, cols):
                M[r][j] -= factor * M[rank][j]

        rank += 1

    return rank


# =============================================================================
# 9. COMPLEXITY
# =============================================================================

def complexity(values):
    if not values:
        return (0, 0, 0)

    primitive = primitive_integer_signature(values)

    max_num_digits = 0
    max_den_digits = 0
    total_digits = 0

    for x in values:
        max_num_digits = max(max_num_digits, len(str(abs(x.numerator))))
        max_den_digits = max(max_den_digits, len(str(x.denominator)))
        total_digits += len(str(abs(x.numerator)))
        total_digits += len(str(x.denominator))

    return (
        max_num_digits,
        max_den_digits,
        total_digits,
        len(primitive),
    )


def primitive_integer_signature(values):
    if not values:
        return []

    common_den = 1

    for x in values:
        common_den = lcm(common_den, x.denominator)

    ints = [
        x.numerator * (common_den // x.denominator)
        for x in values
    ]

    g = 0

    for x in ints:
        g = gcd(g, abs(x))

    if g != 0:
        ints = [x // g for x in ints]

    for x in ints:
        if x != 0:
            if x < 0:
                ints = [-y for y in ints]
            break

    return ints


# =============================================================================
# 10. SMALL-PRIME FACTOR AUDIT
# =============================================================================

def small_prime_factorization(n, bound=1000):
    """
    Exact trial division only up to `bound`.
    Returns:
        factors,
        remaining_cofactor

    A remaining cofactor is NOT claimed prime.
    """
    n = abs(n)

    factors = []

    if n == 0:
        return [], 0

    p = 2

    while p <= bound and p * p <= n:
        exponent = 0

        while n % p == 0:
            n //= p
            exponent += 1

        if exponent:
            factors.append((p, exponent))

        if p == 2:
            p = 3
        else:
            p += 2

    return factors, n


def smoothness_report(n):
    factors, remaining = small_prime_factorization(n)

    return {
        "small_factors": factors,
        "remaining": remaining,
        "completely_factored_to_1000": remaining == 1,
    }


# =============================================================================
# 11. MINIMUM ABSOLUTE DETERMINANT / LOG SCALE
# =============================================================================

def decimal_digits_integer(n):
    n = abs(n)

    if n == 0:
        return 1

    return len(str(n))


def safe_log10_integer(n):
    n = abs(n)

    if n == 0:
        return float("-inf")

    return len(str(n)) - 1


# =============================================================================
# 12. VERIFY OPERATOR POINTWISE
# =============================================================================

def eval_operator(solution, p, r):
    d = D[p] - r
    basis = monomial_basis_value(p, d)

    f0 = sum(
        solution[i] * basis[i]
        for i in range(6)
    )

    f1 = sum(
        solution[6 + i] * basis[i]
        for i in range(6)
    )

    h0 = (
        solution[12]
        + solution[13] * p
    )

    result = (
        f0 * q_value(p, r)
        + f1 * q_value(p, r + 1)
    )

    if d == 0:
        result += h0

    return result


def verify_operator(solution):
    results = []

    for p, p_next in [(1, 3), (3, 5), (5, 7)]:
        for r in range(len(Q[p])):
            lhs = q_value(p_next, r)
            rhs = eval_operator(solution, p, r)

            results.append(
                {
                    "p": p,
                    "p_next": p_next,
                    "r": r,
                    "exact": lhs == rhs,
                }
            )

    return results


# =============================================================================
# 13. CRAMER COLUMN DETERMINANTS
# =============================================================================

def replace_column(A, col, b):
    M = [row[:] for row in A]

    for i in range(len(M)):
        M[i][col] = b[i]

    return M


def compute_cramer_determinants(integer_A, integer_b):
    base_det = bareiss_det(integer_A)

    numerators = []

    for col in range(len(integer_A[0])):
        replaced = replace_column(
            integer_A,
            col,
            integer_b,
        )

        numerators.append(
            bareiss_det(replaced)
        )

    return base_det, numerators


# =============================================================================
# 14. MAIN
# =============================================================================

def main():
    print("=" * 78)
    print("EXPERIMENT 131 — EXACT B-ODD OPERATOR DETERMINANT / CRAMER-SCALE AUDIT")
    print("=" * 78)
    print()

    # -------------------------------------------------------------------------
    # 1. Data validation
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("1. EXACT DATA VALIDATION")
    print("=" * 78)

    data_ok = falling_values()

    print("  p=1: entries=%d degree=%d" % (len(Q[1]), len(Q[1]) - 1))
    print("  p=3: entries=%d degree=%d" % (len(Q[3]), len(Q[3]) - 1))
    print("  p=5: entries=%d degree=%d" % (len(Q[5]), len(Q[5]) - 1))
    print("  p=7: entries=%d degree=%d" % (len(Q[7]), len(Q[7]) - 1))
    print("  data_exact=%s" % data_ok)
    print()

    # -------------------------------------------------------------------------
    # 2. Build system
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("2. EXACT 14 x 14 B-ODD SYSTEM")
    print("=" * 78)

    A, b, meta = build_system()

    print("  rows=%d" % len(A))
    print("  columns=%d" % len(A[0]))
    print("  unknowns=%d" % len(UNKNOWN_NAMES))
    print("  rank=%d" % rank_exact(A))
    print()

    # -------------------------------------------------------------------------
    # 3. Exact solve
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("3. EXACT OPERATOR RECONSTRUCTION")
    print("=" * 78)

    solution = solve_exact(A, b)

    operator_exact = solution is not None

    print("  unique_exact=%s" % operator_exact)

    if not operator_exact:
        print()
        print("OPERATOR RECONSTRUCTION FAILED")
        return

    print()

    for name, value in zip(UNKNOWN_NAMES, solution):
        print("  %-20s = %s" % (name, value))

    print()

    # -------------------------------------------------------------------------
    # 4. Pointwise verification
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("4. POINTWISE RECONSTRUCTION")
    print("=" * 78)

    verification = verify_operator(solution)

    pointwise_ok = True

    for item in verification:
        print(
            "  %d->%d r=%d exact=%s"
            % (
                item["p"],
                item["p_next"],
                item["r"],
                item["exact"],
            )
        )

        pointwise_ok = pointwise_ok and item["exact"]

    print()

    # -------------------------------------------------------------------------
    # 5. Clear denominators row-by-row
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("5. EXACT INTEGER SYSTEM AFTER ROW CLEARING")
    print("=" * 78)

    integer_A, integer_b, row_scale_product = row_clear_fraction_matrix(A, b)

    integer_rank = rank_exact(integer_A)

    print("  integer_rank=%d" % integer_rank)
    print("  row_scale_product_digits=%d" %
          decimal_digits_integer(row_scale_product))

    print()

    # -------------------------------------------------------------------------
    # 6. Determinant
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("6. EXACT DETERMINANT AUDIT")
    print("=" * 78)

    det_A = bareiss_det(integer_A)

    print("  det(integer_A)=...")
    print("  determinant_digits=%d" % decimal_digits_integer(det_A))
    print("  determinant_zero=%s" % (det_A == 0))
    print("  determinant_sign=%s" % ("negative" if det_A < 0 else "positive"))

    print()

    # Exact determinant of the ORIGINAL rational matrix:
    det_fraction = Fraction(det_A, row_scale_product)

    print(
        "  det(A_rational) numerator_digits=%d"
        % decimal_digits_integer(det_fraction.numerator)
    )

    print(
        "  det(A_rational) denominator_digits=%d"
        % decimal_digits_integer(det_fraction.denominator)
    )

    print()

    # -------------------------------------------------------------------------
    # 7. Cramer's rule audit
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("7. EXACT CRAMER DETERMINANT AUDIT")
    print("=" * 78)

    cramer_denominator = det_A

    cramer_numerators = []

    if len(integer_A) == 14 and len(integer_A[0]) == 14:
        cramer_base, cramer_numerators = compute_cramer_determinants(
            integer_A,
            integer_b,
        )

        print(
            "  base_determinant_matches=%s"
            % (cramer_base == det_A)
        )

        for i, numerator in enumerate(cramer_numerators):
            print(
                "  %-20s numerator_digits=%d"
                % (
                    UNKNOWN_NAMES[i],
                    decimal_digits_integer(numerator),
                )
            )

        print()
    else:
        print("  skipped — system is not 14 x 14")

    # -------------------------------------------------------------------------
    # 8. Compare Cramer quotient with solved coefficient
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("8. CRAMER-QUOTIENT EXACTNESS")
    print("=" * 78)

    cramer_exact = True

    for i, numerator in enumerate(cramer_numerators):
        candidate = Fraction(numerator, cramer_denominator)

        ok = candidate == solution[i]

        print(
            "  %-20s exact=%s"
            % (
                UNKNOWN_NAMES[i],
                ok,
            )
        )

        cramer_exact = cramer_exact and ok

    print()

    # -------------------------------------------------------------------------
    # 9. Denominator comparison
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("9. COEFFICIENT-DENOMINATOR VS DETERMINANT SCALE")
    print("=" * 78)

    all_denominators = [
        abs(x.denominator)
        for x in solution
        if x != 0
    ]

    max_coeff_den = max(all_denominators)

    print(
        "  max_coefficient_denominator_digits=%d"
        % decimal_digits_integer(max_coeff_den)
    )

    print(
        "  integer_system_determinant_digits=%d"
        % decimal_digits_integer(abs(det_A))
    )

    print(
        "  rational_system_determinant_denominator_digits=%d"
        % decimal_digits_integer(det_fraction.denominator)
    )

    print()

    # Does determinant denominator absorb every coefficient denominator?
    denominator_absorption = True

    for denom in all_denominators:
        if det_A % denom != 0:
            denominator_absorption = False
            break

    print(
        "  |det(integer_A)| divisible_by_all_coefficient_denominators=%s"
        % denominator_absorption
    )

    print()

    # -------------------------------------------------------------------------
    # 10. Common coefficient denominator
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("10. COMMON COEFFICIENT DENOMINATOR")
    print("=" * 78)

    common_den = 1

    for x in solution:
        common_den = lcm(common_den, x.denominator)

    print(
        "  common_denominator_digits=%d"
        % decimal_digits_integer(common_den)
    )

    print("  common_denominator=%s" % common_den)

    print(
        "  determinant_divisible_by_common_denominator=%s"
        % (det_A % common_den == 0)
    )

    print()

    # -------------------------------------------------------------------------
    # 11. Small-prime smoothness
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("11. SMALL-PRIME SMOOTHNESS AUDIT")
    print("=" * 78)

    smooth_det = smoothness_report(det_A)

    print(
        "  determinant_small_prime_factors=%s"
        % smooth_det["small_factors"]
    )

    print(
        "  determinant_remaining_after_primes_le_1000=%s"
        % smooth_det["remaining"]
    )

    print(
        "  completely_factored_to_1000=%s"
        % smooth_det["completely_factored_to_1000"]
    )

    smooth_den = smoothness_report(common_den)

    print(
        "  coefficient_denominator_small_prime_factors=%s"
        % smooth_den["small_factors"]
    )

    print(
        "  coefficient_denominator_remaining_after_primes_le_1000=%s"
        % smooth_den["remaining"]
    )

    print(
        "  coefficient_denominator_completely_factored_to_1000=%s"
        % smooth_den["completely_factored_to_1000"]
    )

    print()

    # -------------------------------------------------------------------------
    # 12. Primitive signatures
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("12. PRIMITIVE OPERATOR SIGNATURES")
    print("=" * 78)

    blocks = {
        "F0": solution[0:6],
        "F1": solution[6:12],
        "H0": solution[12:14],
    }

    for name, coeffs in blocks.items():
        signature = primitive_integer_signature(coeffs)

        print("  %s:" % name)
        print("    primitive_signature=%s" % signature)
        print(
            "    complexity=%s"
            % (complexity(coeffs),)
        )

    print()

    # -------------------------------------------------------------------------
    # 13. Fact-check friend's determinant hypothesis
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("13. DETERMINANT-HYPOTHESIS FACT CHECK")
    print("=" * 78)

    print(
        """
The exact 14-unknown system is genuinely full rank, so Cramer's-rule
denominators are mathematically relevant.

However, the following stronger claims must be checked separately:

  (a) Are the ~70-digit integers actually determinant/cofactor scales?
  (b) Does the determinant dominate the coefficient denominators?
  (c) Are the determinants/factors unusually smooth?
  (d) Does a basis change dramatically reduce the rational complexity?

This experiment records the exact quantities rather than assuming the answer.
"""
    )

    determinant_scale_explains = (
        decimal_digits_integer(abs(det_A))
        >= decimal_digits_integer(max_coeff_den)
    )

    print(
        "  determinant_scale_at_least_coefficient_denominator_scale=%s"
        % determinant_scale_explains
    )

    print(
        "  cramer_rule_exact=%s"
        % cramer_exact
    )

    print()

    # -------------------------------------------------------------------------
    # 14. Final exactness
    # -------------------------------------------------------------------------
    print("=" * 78)
    print("14. FINAL EXACTNESS")
    print("=" * 78)

    failures = 0

    if not data_ok:
        failures += 1

    if not operator_exact:
        failures += 1

    if rank_exact(A) != 14:
        failures += 1

    if not pointwise_ok:
        failures += 1

    if not cramer_exact:
        failures += 1

    print("  data_exact=%s" % data_ok)
    print("  operator_exact=%s" % operator_exact)
    print("  rank_14_exact=%s" % (rank_exact(A) == 14))
    print("  pointwise_exact=%s" % pointwise_ok)
    print("  cramer_exact=%s" % cramer_exact)
    print("  failures=%d" % failures)

    if failures == 0:
        print("  ALL BASIC CHECKS PASS=True")
    else:
        print("  ALL BASIC CHECKS PASS=False")

    print()
    print("EXPERIMENT 131 COMPLETE")


if __name__ == "__main__":
    main()

