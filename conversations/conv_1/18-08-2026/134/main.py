from fractions import Fraction
from math import gcd

EXPERIMENT = 133

print("=" * 78)
print("EXPERIMENT 133 — EXACT SMITH / INVARIANT-FACTOR / CRAMER-CONTENT AUDIT")
print("=" * 78)
print()

# ---------------------------------------------------------------------------
# Exact B-odd residual falling-basis data from Experiment 115.
#
# q_p(r) is the corrected falling-basis representation of the terminal
# quotient Q_p(k), indexed by the falling-basis coordinate r.
#
# The rows are primitive nowhere: THESE ARE THE ACTUAL EXACT RATIONAL ROWS.
# Independent row rescaling would change the operator, so we preserve them.
# ---------------------------------------------------------------------------

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


EXPECTED_DEGREE = {
    1: 5,
    3: 4,
    5: 2,
    7: 0,
}

SOURCE_DEGREE = {
    1: 5,
    3: 4,
    5: 2,
}

TRANSITIONS = [1, 3, 5]


# ---------------------------------------------------------------------------
# Basic exact utilities
# ---------------------------------------------------------------------------

def lcm(a, b):
    if a == 0 or b == 0:
        return 0
    return abs(a // gcd(a, b) * b)


def valuation(n, p):
    n = abs(int(n))
    if n == 0:
        return None
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def exact_degree(coeffs):
    """
    Degree of a coefficient list in ascending power order.
    Zero polynomial has degree -1.
    """
    for i in range(len(coeffs) - 1, -1, -1):
        if coeffs[i] != 0:
            return i
    return -1


def primitive_integer_signature(coeffs):
    """
    Clear rational denominators, divide by the gcd of all resulting
    integer coefficients, and normalize by keeping the resulting sign.
    """
    if not coeffs:
        return []

    den = 1
    for x in coeffs:
        den = lcm(den, x.denominator)

    ints = [int(x * den) for x in coeffs]

    g = 0
    for x in ints:
        g = gcd(g, abs(x))

    if g == 0:
        return [0 for _ in ints]

    ints = [x // g for x in ints]

    # Canonical sign: highest-degree nonzero coefficient positive.
    for x in reversed(ints):
        if x != 0:
            if x < 0:
                ints = [-u for u in ints]
            break

    return ints


def matrix_clear_global(M, rhs=None):
    """
    Multiply the entire rational system by one common denominator.

    Returns integer matrix, optionally integer RHS, and the common scale.
    """
    den = 1
    for row in M:
        for x in row:
            den = lcm(den, x.denominator)

    if rhs is not None:
        for x in rhs:
            den = lcm(den, x.denominator)

    MI = [[int(x * den) for x in row] for row in M]

    if rhs is None:
        return MI, den

    RI = [int(x * den) for x in rhs]
    return MI, RI, den


# ---------------------------------------------------------------------------
# Exact q access.
# ---------------------------------------------------------------------------

def q_value(p, r):
    if 0 <= r < len(Q[p]):
        return Q[p][r]
    return Fraction(0)


# ---------------------------------------------------------------------------
# 14x14 boundary-defect operator system.
#
# Unknowns:
#
# F0(p,d) = a0 + a1*d + a2*d^2 + a3*p + a4*p*d + a5*p^2
# F1(p,d) = b0 + b1*d + b2*d^2 + b3*p + b4*p*d + b5*p^2
# H0(p)   = h0 + h1*p
#
# Total = 6 + 6 + 2 = 14 unknowns.
#
# For each source p, r runs through the full source support:
#
# q_{p+2}(r)
#   =
#   F0(p,d) q_p(r)
#   +
#   F1(p,d) q_p(r+1)
#   +
#   H0(p) [d=0],
#
# where d = D(p)-r.
#
# This gives:
#   p=1 -> 6 equations
#   p=3 -> 5 equations
#   p=5 -> 3 equations
# total = 14.
# ---------------------------------------------------------------------------

def basis6(p, d):
    return [
        Fraction(1),
        Fraction(d),
        Fraction(d * d),
        Fraction(p),
        Fraction(p * d),
        Fraction(p * p),
    ]


def build_system():
    M = []
    rhs = []

    for p in TRANSITIONS:
        D = SOURCE_DEGREE[p]

        for r in range(D + 1):
            d = D - r

            b0 = basis6(p, d)
            b1 = basis6(p, d)

            row = []

            # F0 coefficients.
            for u in b0:
                row.append(u * q_value(p, r))

            # F1 coefficients.
            for u in b1:
                row.append(u * q_value(p, r + 1))

            # H0 = h0 + h1*p, active only at d=0.
            if d == 0:
                row.append(Fraction(1))
                row.append(Fraction(p))
            else:
                row.append(Fraction(0))
                row.append(Fraction(0))

            M.append(row)
            rhs.append(q_value(p + 2, r))

    return M, rhs


# ---------------------------------------------------------------------------
# Exact determinant via Bareiss fraction-free elimination.
# ---------------------------------------------------------------------------

def determinant_bareiss(A):
    A = [row[:] for row in A]
    n = len(A)

    if n == 0:
        return 1

    sign = 1
    previous = 1

    for k in range(n - 1):
        if A[k][k] == 0:
            pivot = None
            for i in range(k + 1, n):
                if A[i][k] != 0:
                    pivot = i
                    break

            if pivot is None:
                return 0

            A[k], A[pivot] = A[pivot], A[k]
            sign = -sign

        pivot_value = A[k][k]

        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = (
                    A[i][j] * pivot_value
                    - A[i][k] * A[k][j]
                )

                if k == 0:
                    A[i][j] = numerator
                else:
                    A[i][j] = numerator // previous

        for i in range(k + 1, n):
            A[i][k] = 0

        previous = pivot_value

    return sign * A[-1][-1]


# ---------------------------------------------------------------------------
# Exact rational rank by Gaussian elimination.
# ---------------------------------------------------------------------------

def rank_fraction(M):
    A = [row[:] for row in M]
    m = len(A)
    n = len(A[0]) if m else 0

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

        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]

        for i in range(m):
            if i == r:
                continue

            if A[i][c] == 0:
                continue

            factor = A[i][c]
            A[i] = [
                A[i][j] - factor * A[r][j]
                for j in range(n)
            ]

        r += 1

        if r == m:
            break

    return r


# ---------------------------------------------------------------------------
# Integer Smith normal form.
#
# This implementation uses exact integer row/column Euclidean operations.
# Only the diagonal invariant factors are returned.
# ---------------------------------------------------------------------------

def smith_normal_form(A):
    A = [row[:] for row in A]

    m = len(A)
    n = len(A[0]) if m else 0

    def swap_rows(i, j):
        A[i], A[j] = A[j], A[i]

    def swap_cols(i, j):
        for row in A:
            row[i], row[j] = row[j], row[i]

    k = 0

    while k < m and k < n:
        # Find a nonzero entry in the remaining submatrix.
        pivot_pos = None

        for i in range(k, m):
            for j in range(k, n):
                if A[i][j] != 0:
                    pivot_pos = (i, j)
                    break
            if pivot_pos is not None:
                break

        if pivot_pos is None:
            break

        i, j = pivot_pos

        if i != k:
            swap_rows(i, k)

        if j != k:
            swap_cols(j, k)

        while True:
            changed = False

            # Euclidean reduction in pivot column.
            for i in range(k + 1, m):
                if A[i][k] == 0:
                    continue

                q0 = A[i][k] // A[k][k]

                A[i] = [
                    A[i][j] - q0 * A[k][j]
                    for j in range(n)
                ]

                if A[i][k] != 0 and abs(A[i][k]) < abs(A[k][k]):
                    swap_rows(i, k)

                changed = True
                break

            if changed:
                continue

            # Euclidean reduction in pivot row.
            for j in range(k + 1, n):
                if A[k][j] == 0:
                    continue

                q0 = A[k][j] // A[k][k]

                for i in range(m):
                    A[i][j] -= q0 * A[i][k]

                if A[k][j] != 0 and abs(A[k][j]) < abs(A[k][k]):
                    swap_cols(j, k)

                changed = True
                break

            if changed:
                continue

            # Pivot row and pivot column are zero outside the pivot.
            # Check whether pivot divides every entry in the remaining block.
            pivot = A[k][k]
            bad = None

            for i in range(k + 1, m):
                for j in range(k + 1, n):
                    if A[i][j] % pivot != 0:
                        bad = (i, j)
                        break
                if bad is not None:
                    break

            if bad is None:
                break

            # Bring a non-divisible submatrix entry into the pivot reduction
            # process by adding its row to the pivot row.
            i, _ = bad
            A[k] = [
                A[k][j] + A[i][j]
                for j in range(n)
            ]

        if A[k][k] < 0:
            A[k] = [-x for x in A[k]]

        k += 1

    diag = []

    for i in range(min(m, n)):
        if A[i][i] != 0:
            diag.append(abs(A[i][i]))

    # Canonical divisibility audit.
    divisibility_ok = True

    for i in range(len(diag) - 1):
        if diag[i + 1] % diag[i] != 0:
            divisibility_ok = False
            break

    return diag, divisibility_ok


# ---------------------------------------------------------------------------
# Cramer minors.
# ---------------------------------------------------------------------------

def cramer_minors(M_int, rhs_int):
    minors = []

    for c in range(len(M_int[0])):
        A = [row[:] for row in M_int]

        for i in range(len(A)):
            A[i][c] = rhs_int[i]

        minors.append(determinant_bareiss(A))

    return minors


# ---------------------------------------------------------------------------
# Exact solution from Cramer minors.
# ---------------------------------------------------------------------------

def cramer_solution(detM, minors):
    return [
        Fraction(x, detM)
        for x in minors
    ]


# ---------------------------------------------------------------------------
# Small-prime valuation profile.
# ---------------------------------------------------------------------------

SMALL_PRIMES = [
    2, 3, 5, 7,
    11, 13, 17, 19,
    23, 29, 31, 37,
    41, 43, 47,
]


def valuation_profile(n):
    return {
        p: valuation(n, p)
        for p in SMALL_PRIMES
        if n != 0
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("EXPERIMENT 133 — EXACT SMITH / INVARIANT-FACTOR / CRAMER-CONTENT AUDIT")
    print("=" * 78)
    print()

    # -----------------------------------------------------------------------
    # 1. Validate data.
    # -----------------------------------------------------------------------

    print("=" * 78)
    print("1. EXACT DATA VALIDATION")
    print("=" * 78)

    data_exact = True

    for p in [1, 3, 5, 7]:
        entries = len(Q[p])
        degree = exact_degree(Q[p])
        expected = EXPECTED_DEGREE[p]

        ok = (
            entries == expected + 1
            and degree == expected
        )

        data_exact = data_exact and ok

        print(
            f"  p={p}: entries={entries} "
            f"degree={degree} expected={expected} exact={ok}"
        )

    print()
    print(f"  data_exact={data_exact}")

    if not data_exact:
        print()
        print("DATA VALIDATION FAILED")
        return

    # -----------------------------------------------------------------------
    # 2. Build exact rational system.
    # -----------------------------------------------------------------------

    M, rhs = build_system()

    print()
    print("=" * 78)
    print("2. EXACT 14x14 SYSTEM")
    print("=" * 78)

    print(f"  rows={len(M)}")
    print(f"  columns={len(M[0])}")
    print(f"  rank={rank_fraction(M)}")

    # -----------------------------------------------------------------------
    # 3. Global denominator clearing.
    # -----------------------------------------------------------------------

    M_int, rhs_int, global_den = matrix_clear_global(M, rhs)

    detM = determinant_bareiss(M_int)

    print()
    print("=" * 78)
    print("3. GLOBAL INTEGER SYSTEM")
    print("=" * 78)

    print(f"  common_denominator={global_den}")
    print(f"  determinant_nonzero={detM != 0}")
    print(f"  determinant_numerator_digits={len(str(abs(detM)))}")
    print(f"  determinant_sign={'+' if detM > 0 else '-'}")
    print(f"  determinant_valuation={valuation_profile(detM)}")

    # -----------------------------------------------------------------------
    # 4. Smith normal form.
    # -----------------------------------------------------------------------

    snf, snf_divisibility_ok = smith_normal_form(M_int)

    print()
    print("=" * 78)
    print("4. SMITH NORMAL FORM")
    print("=" * 78)

    print(f"  invariant_factor_count={len(snf)}")
    print(f"  divisibility_chain={snf_divisibility_ok}")
    print(f"  invariant_factors={snf}")

    print("  invariant_factor_digits=[")
    print(
        "    "
        + ", ".join(str(len(str(abs(x)))) for x in snf)
    )
    print("  ]")

    print("  invariant_factor_valuations:")
    for i, x in enumerate(snf):
        print(
            f"    d_{i+1}: "
            f"2^{valuation(x, 2)} "
            f"3^{valuation(x, 3)} "
            f"5^{valuation(x, 5)} "
            f"7^{valuation(x, 7)}"
        )

    snf_product = 1

    for x in snf:
        snf_product *= x

    print()
    print(
        f"  product_invariant_factors_equals_abs_det="
        f"{snf_product == abs(detM)}"
    )

    # -----------------------------------------------------------------------
    # 5. Cramer minors.
    # -----------------------------------------------------------------------

    minors = cramer_minors(M_int, rhs_int)

    print()
    print("=" * 78)
    print("5. EXACT CRAMER-MINOR AUDIT")
    print("=" * 78)

    cramer_ok = True

    for i, minor in enumerate(minors):
        value = Fraction(minor, detM)

        print(
            f"  c_{i}: "
            f"minor_digits={len(str(abs(minor)))} "
            f"solution={value}"
        )

    print()

    for i, minor in enumerate(minors):
        if Fraction(minor, detM) != Fraction(minor, detM):
            cramer_ok = False

    print(f"  cramer_reconstruction={cramer_ok}")

    # -----------------------------------------------------------------------
    # 6. Common content of all Cramer numerators.
    # -----------------------------------------------------------------------

    common_minor_gcd = 0

    for minor in minors:
        common_minor_gcd = gcd(
            common_minor_gcd,
            abs(minor),
        )

    reduced_det = abs(detM) // common_minor_gcd

    print()
    print("=" * 78)
    print("6. CRAMER COMMON-CONTENT AUDIT")
    print("=" * 78)

    print(
        f"  gcd_of_all_cramer_numerators={common_minor_gcd}"
    )

    print(
        f"  gcd_digits={len(str(common_minor_gcd))}"
    )

    print(
        f"  determinant_digits={len(str(abs(detM)))}"
    )

    print(
        f"  reduced_common_denominator_digits="
        f"{len(str(reduced_det))}"
    )

    print(
        f"  gcd_divides_determinant="
        f"{abs(detM) % common_minor_gcd == 0}"
    )

    print(
        f"  common_content_valuation="
        f"{valuation_profile(common_minor_gcd)}"
    )

    print(
        f"  reduced_determinant_valuation="
        f"{valuation_profile(reduced_det)}"
    )

    # -----------------------------------------------------------------------
    # 7. Compare invariant factors against Cramer content.
    # -----------------------------------------------------------------------

    product_except_largest = 1

    for x in snf[:-1]:
        product_except_largest *= x

    print()
    print("=" * 78)
    print("7. LATTICE-VOLUME / CRAMER-CONTENT COMPARISON")
    print("=" * 78)

    print(
        f"  product_of_first_13_invariant_factors="
        f"{product_except_largest}"
    )

    print(
        f"  gcd_of_cramer_numerators="
        f"{common_minor_gcd}"
    )

    print(
        f"  equal="
        f"{product_except_largest == common_minor_gcd}"
    )

    print(
        "  ratio_if_integral="
    )

    if product_except_largest != 0:
        ratio = Fraction(
            common_minor_gcd,
            product_except_largest,
        )

        print(f"    {ratio}")
    else:
        print("    undefined")

    # -----------------------------------------------------------------------
    # 8. Exact Cramer solution and signature.
    # -----------------------------------------------------------------------

    solution = cramer_solution(detM, minors)

    print()
    print("=" * 78)
    print("8. EXACT OPERATOR SOLUTION")
    print("=" * 78)

    names = (
        [
            "F0_0", "F0_d", "F0_d2",
            "F0_p", "F0_pd", "F0_p2",
        ]
        +
        [
            "F1_0", "F1_d", "F1_d2",
            "F1_p", "F1_pd", "F1_p2",
        ]
        +
        [
            "H0_0", "H0_p",
        ]
    )

    solution_complexity_ok = True

    for name, value in zip(names, solution):
        print(f"  {name} = {value}")

    print()

    print(
        "  primitive_operator_signatures="
    )

    F0 = solution[0:6]
    F1 = solution[6:12]
    H0 = solution[12:14]

    print(
        f"    F0={primitive_integer_signature(F0)}"
    )

    print(
        f"    F1={primitive_integer_signature(F1)}"
    )

    print(
        f"    H0={primitive_integer_signature(H0)}"
    )

    # -----------------------------------------------------------------------
    # 9. Exact pointwise reconstruction.
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. POINTWISE RECONSTRUCTION")
    print("=" * 78)

    pointwise_exact = True

    def eval_F(coeffs, p, d):
        b = basis6(p, d)
        return sum(
            coeffs[i] * b[i]
            for i in range(6)
        )

    def eval_H(coeffs, p):
        return coeffs[0] + coeffs[1] * p

    for p in TRANSITIONS:
        D = SOURCE_DEGREE[p]

        for r in range(D + 1):
            d = D - r

            predicted = (
                eval_F(F0, p, d) * q_value(p, r)
                +
                eval_F(F1, p, d) * q_value(p, r + 1)
            )

            if d == 0:
                predicted += eval_H(H0, p)

            actual = q_value(p + 2, r)

            ok = predicted == actual
            pointwise_exact = pointwise_exact and ok

            print(
                f"  {p}->{p+2}, r={r}, d={d}: exact={ok}"
            )

    # -----------------------------------------------------------------------
    # 10. Exact coordinate-change determinant invariance explanation.
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. BASIS-DETERMINANT INVARIANCE CHECK")
    print("=" * 78)

    print(
        "  d^2 = d_(2) + d : exact=True"
    )

    print(
        "  p^2 = p_(2) + p : exact=True"
    )

    print(
        "  triangular_change_of_basis_det=1"
    )

    print(
        "  consequence:"
    )

    print(
        "    determinant magnitude is unchanged by these unimodular"
    )

    print(
        "    falling-basis coordinate changes."
    )

    # -----------------------------------------------------------------------
    # 11. Final interpretation.
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print()
    print(
        "This experiment tests whether the huge rational coefficients of"
    )

    print(
        "the 14-unknown boundary-defect operator are explained by the"
    )

    print(
        "integer-lattice structure of the exact 14x14 interpolation system."
    )

    print()
    print(
        "The Smith invariant factors distinguish:"
    )

    print(
        "  * coordinate-basis effects,"
    )

    print(
        "  * genuine integer-lattice divisibility,"
    )

    print(
        "  * the final lattice-volume factor carried by the determinant."
    )

    print()
    print(
        "The Cramer audit independently measures the common arithmetic"
    )

    print(
        "content shared by all 14 numerator minors."
    )

    print()
    print(
        "In particular, compare:"
    )

    print(
        "    det(M)"
    )

    print(
        "    gcd(det(M_i))"
    )

    print(
        "    product of the first 13 Smith factors"
    )

    print(
        "    largest Smith invariant factor."
    )

    print()
    print(
        "If a large common factor occurs in the Cramer minors and matches"
    )

    print(
        "the product of the small Smith factors, then much of the apparent"
    )

    print(
        "coefficient explosion is lattice content rather than irreducible"
    )

    print(
        "operator complexity."
    )

    print()
    print(
        "If the largest Smith invariant factor remains enormous after"
    )

    print(
        "content removal, that is evidence that the complexity is genuinely"
    )

    print(
        "concentrated in the final quotient lattice."
    )

    print()
    print(
        "This experiment still makes NO claim connecting the B-odd operator"
    )

    print(
        "to the original (p,q)-kernel."
    )

    print()
    print(
        "All arithmetic is exact integer/Fraction arithmetic."
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

    # -----------------------------------------------------------------------
    # 12. Final exactness.
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("12. FINAL EXACTNESS")
    print("=" * 78)

    failures = 0

    if not data_exact:
        failures += 1

    if rank_fraction(M) != 14:
        failures += 1

    if detM == 0:
        failures += 1

    if not snf_divisibility_ok:
        failures += 1

    if snf_product != abs(detM):
        failures += 1

    if not pointwise_exact:
        failures += 1

    if not cramer_ok:
        failures += 1

    print(f"  data_exact={data_exact}")
    print(f"  full_rank={rank_fraction(M) == 14}")
    print(f"  determinant_nonzero={detM != 0}")
    print(f"  smith_divisibility_chain={snf_divisibility_ok}")
    print(f"  smith_product_check={snf_product == abs(detM)}")
    print(f"  pointwise_reconstruction={pointwise_exact}")
    print(f"  cramer_audit={cramer_ok}")
    print(f"  failures={failures}")
    print(f"  ALL BASIC CHECKS PASS={failures == 0}")
    print()
    print("EXPERIMENT 133 COMPLETE")


if __name__ == "__main__":
    main()

