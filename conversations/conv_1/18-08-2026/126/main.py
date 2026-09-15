#!/usr/bin/env python3

from fractions import Fraction
from math import gcd


# ============================================================================
# EXPERIMENT 126
# EXACT BOUNDARY-DEFECT FALLING-BASIS / COEFFICIENT NORMALIZATION AUDIT
# ============================================================================

B_ODD = {
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

D = {p: len(row) - 1 for p, row in B_ODD.items()}


# ============================================================================
# BASIC POLYNOMIAL UTILITIES
# ============================================================================

def trim(poly):
    poly = list(poly)
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def poly_degree(poly):
    poly = trim(poly)
    if len(poly) == 1 and poly[0] == 0:
        return -1
    return len(poly) - 1


def poly_eval(poly, x):
    total = Fraction(0)
    for c in reversed(poly):
        total = total * x + c
    return total


# ============================================================================
# EXACT LINEAR ALGEBRA
# ============================================================================

def rref(A):
    if not A:
        return [], []

    M = [[Fraction(x) for x in row] for row in A]
    rows = len(M)
    cols = len(M[0])

    rank = 0
    pivots = []

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

        pivots.append(col)
        rank += 1

        if rank == rows:
            break

    return M, pivots


def matrix_rank(A):
    if not A:
        return 0
    _, pivots = rref(A)
    return len(pivots)


def solve_unique(A, b):
    if not A:
        return None

    nvars = len(A[0])

    aug = [
        [Fraction(x) for x in A[i]] + [Fraction(b[i])]
        for i in range(len(A))
    ]

    R, pivots = rref(aug)

    # Inconsistency
    for row in R:
        if all(row[j] == 0 for j in range(nvars)):
            if row[nvars] != 0:
                return None

    # Need unique solution
    if len(pivots) != nvars:
        return None

    sol = [Fraction(0)] * nvars

    for i, pivot in enumerate(pivots):
        sol[pivot] = R[i][nvars]

    return sol


# ============================================================================
# FALLING FACTORIALS
# ============================================================================

def falling_poly(n):
    if n == 0:
        return [Fraction(1)]

    out = [Fraction(1)]

    for t in range(n):
        nxt = [Fraction(0)] * (len(out) + 1)

        for i, c in enumerate(out):
            nxt[i] -= t * c
            nxt[i + 1] += c

        out = trim(nxt)

    return out


# ============================================================================
# DATA ACCESS
# ============================================================================

def get_q(p, r):
    row = B_ODD[p]

    if 0 <= r < len(row):
        return row[r]

    return Fraction(0)


# ============================================================================
# OPERATOR BASIS
# ============================================================================

def pd_basis(p, d):
    # 1, d, d^2, p, p*d, p^2
    return [
        Fraction(1),
        Fraction(d),
        Fraction(d * d),
        Fraction(p),
        Fraction(p * d),
        Fraction(p * p),
    ]


def p_basis(p):
    return [
        Fraction(1),
        Fraction(p),
    ]


# ============================================================================
# EXPERIMENT-125 OPERATOR RECONSTRUCTION
# ============================================================================

def build_operator_system():
    """
    Model:

      q_{p+2}(r)
        =
        F0(p,d) q_p(r)
        + F1(p,d) q_p(r+1)
        + H0(p) [d=0]

    d = D(p)-r

    F0, F1: total degree <= 2 in (p,d)
    H0: degree <= 1 in p
    """

    A = []
    b = []

    for p in sorted(B_ODD):
        p2 = p + 2

        if p2 not in B_ODD:
            continue

        Dp = D[p]

        # Include several points beyond the visible support.
        # Missing q-values are exactly zero.
        max_r = max(D[p], D[p2]) + 5

        for r in range(max_r + 1):

            d = Dp - r

            basis = pd_basis(p, d)
            pb = p_basis(p)

            equation = [Fraction(0)] * 14

            q0 = get_q(p, r)
            if q0 != 0:
                for j in range(6):
                    equation[j] += q0 * basis[j]

            q1 = get_q(p, r + 1)
            if q1 != 0:
                for j in range(6):
                    equation[6 + j] += q1 * basis[j]

            if d == 0:
                for j in range(2):
                    equation[12 + j] += pb[j]

            A.append(equation)
            b.append(get_q(p2, r))

    return A, b


def reconstruct_operator():
    A, b = build_operator_system()

    rank_A = matrix_rank(A)
    rank_aug = matrix_rank(
        [
            A[i] + [b[i]]
            for i in range(len(A))
        ]
    )

    if rank_A != rank_aug:
        return None, rank_A, rank_aug

    solution = solve_unique(A, b)

    return solution, rank_A, rank_aug


# ============================================================================
# OPERATOR EVALUATION
# ============================================================================

def evaluate_operator(solution, p, d, r):
    F0 = solution[0:6]
    F1 = solution[6:12]
    H0 = solution[12:14]

    basis = pd_basis(p, d)
    pb = p_basis(p)

    total = Fraction(0)

    total += get_q(p, r) * sum(
        F0[j] * basis[j]
        for j in range(6)
    )

    total += get_q(p, r + 1) * sum(
        F1[j] * basis[j]
        for j in range(6)
    )

    if d == 0:
        total += sum(
            H0[j] * pb[j]
            for j in range(2)
        )

    return total


def pointwise_failures(solution):
    failures = []

    for p in sorted(B_ODD):
        p2 = p + 2

        if p2 not in B_ODD:
            continue

        max_r = max(D[p], D[p2]) + 5

        for r in range(max_r + 1):
            lhs = get_q(p2, r)
            rhs = evaluate_operator(
                solution,
                p,
                D[p] - r,
                r,
            )

            if lhs != rhs:
                failures.append(
                    (p, p2, r, lhs, rhs)
                )

    return failures


# ============================================================================
# COMPLEXITY
# ============================================================================

def primitive_signature(coeffs):
    if not coeffs:
        return []

    den_lcm = 1

    for c in coeffs:
        den_lcm = (
            den_lcm
            * c.denominator
            // gcd(den_lcm, c.denominator)
        )

    integers = [
        int(c * den_lcm)
        for c in coeffs
    ]

    g = 0
    for n in integers:
        g = gcd(g, abs(n))

    if g == 0:
        return [0] * len(integers)

    integers = [n // g for n in integers]

    # Canonical sign: first nonzero entry positive.
    for n in integers:
        if n != 0:
            if n < 0:
                integers = [-x for x in integers]
            break

    return integers


def complexity(coeffs):
    if not coeffs:
        return (0, 0, 0)

    sig = primitive_signature(coeffs)

    max_num_digits = max(
        len(str(abs(c.numerator)))
        for c in coeffs
    )

    max_den_digits = max(
        len(str(c.denominator))
        for c in coeffs
    )

    max_signature_digits = max(
        (
            len(str(abs(int(x))))
            for x in sig
            if x != 0
        ),
        default=0,
    )

    return (
        max_num_digits,
        max_den_digits,
        max_signature_digits,
    )


# ============================================================================
# BIVARIATE MONOMIAL -> FALLING BASIS
# ============================================================================

def monomial_to_falling_2d(block):
    """
    Input ordering:

      [a00, a01, a02, a10, a11, a20]

    representing

      a00
      + a01*d
      + a02*d^2
      + a10*p
      + a11*p*d
      + a20*p^2

    Since

      d^2 = d_(2) + d
      p^2 = p_(2) + p,

    the falling-basis coefficients are

      [a00,
       a01+a02,
       a02,
       a10+a20,
       a11,
       a20].
    """

    a00, a01, a02, a10, a11, a20 = block

    return [
        a00,
        a01 + a02,
        a02,
        a10 + a20,
        a11,
        a20,
    ]


def falling_to_monomial_2d(block):
    b00, b01, b02, b10, b11, b20 = block

    return [
        b00,
        b01 - b02,
        b02,
        b10 - b20,
        b11,
        b20,
    ]


def falling_roundtrip(block):
    return (
        falling_to_monomial_2d(
            monomial_to_falling_2d(block)
        )
        == list(block)
    )


# ============================================================================
# ROOT / PROPORTIONALITY AUDITS
# ============================================================================

def integer_slice_roots(block, bound=20):
    roots = []

    for p in range(-bound, bound + 1):
        value = sum(
            block[j] * pd_basis(p, 0)[j]
            for j in range(6)
        )

        if value == 0:
            roots.append(p)

    return roots


def proportional(block_a, block_b):
    scalar = None

    for a, b in zip(block_a, block_b):

        if b == 0:
            if a != 0:
                return False, None
            continue

        current = a / b

        if scalar is None:
            scalar = current
        elif current != scalar:
            return False, None

    if scalar is None:
        return True, Fraction(0)

    return True, scalar


def boundary_source_roots(H):
    roots = []

    for p in range(-20, 21):
        value = H[0] + H[1] * p

        if value == 0:
            roots.append(p)

    return roots


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 126 — EXACT BOUNDARY-DEFECT "
        "FALLING-BASIS / COEFFICIENT NORMALIZATION AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------------
    # 1
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. DATA VALIDATION")
    print("=" * 78)

    data_exact = True

    for p in sorted(B_ODD):

        ok = (
            isinstance(B_ODD[p], list)
            and all(isinstance(x, Fraction) for x in B_ODD[p])
            and len(B_ODD[p]) == D[p] + 1
        )

        if not ok:
            data_exact = False

        print(
            f"  p={p}: D={D[p]} "
            f"entries={len(B_ODD[p])} "
            f"exact={ok}"
        )

    # ------------------------------------------------------------------------
    # 2
    # ------------------------------------------------------------------------

    solution, rank_A, rank_aug = reconstruct_operator()

    operator_exact = solution is not None

    print()
    print("=" * 78)
    print("2. EXACT OPERATOR RECONSTRUCTION")
    print("=" * 78)

    print(f"  rank={rank_A}")
    print(f"  augmented_rank={rank_aug}")
    print(f"  unknowns=14")
    print(f"  unique_exact={operator_exact}")

    if not operator_exact:

        print()
        print("  OPERATOR RECONSTRUCTION FAILED")

        print()
        print("=" * 78)
        print("FINAL EXACTNESS")
        print("=" * 78)

        overall = (
            data_exact
            and operator_exact
        )

        print(f"  data_exact={data_exact}")
        print(f"  operator_exact={operator_exact}")
        print(f"  ALL BASIC CHECKS PASS={overall}")
        return

    # ------------------------------------------------------------------------
    # 3
    # ------------------------------------------------------------------------

    failures = pointwise_failures(solution)
    pointwise_exact = len(failures) == 0

    F0 = solution[0:6]
    F1 = solution[6:12]
    H0 = solution[12:14]

    print()
    print("=" * 78)
    print("3. RAW OPERATOR COMPLEXITY")
    print("=" * 78)

    for name, block in (
        ("F_0", F0),
        ("F_1", F1),
        ("H_0", H0),
    ):

        print(
            f"  {name}: degree={len(block)-1} "
            f"nonzero={sum(c != 0 for c in block)}"
        )

        print(
            f"    complexity={complexity(block)}"
        )

        print(
            f"    primitive_signature="
            f"{primitive_signature(block)}"
        )

    # ------------------------------------------------------------------------
    # 4
    # ------------------------------------------------------------------------

    fall_F0 = monomial_to_falling_2d(F0)
    fall_F1 = monomial_to_falling_2d(F1)

    H0_falling = list(H0)

    fall_F0_roundtrip = falling_roundtrip(F0)
    fall_F1_roundtrip = falling_roundtrip(F1)

    print()
    print("=" * 78)
    print("4. FALLING-BASIS CONVERSION")
    print("=" * 78)

    print(f"  F_0 falling={fall_F0}")
    print(
        f"    complexity={complexity(fall_F0)}"
    )
    print(
        f"    roundtrip={fall_F0_roundtrip}"
    )

    print(f"  F_1 falling={fall_F1}")
    print(
        f"    complexity={complexity(fall_F1)}"
    )
    print(
        f"    roundtrip={fall_F1_roundtrip}"
    )

    print(f"  H_0 falling={H0_falling}")
    print(
        f"    complexity={complexity(H0_falling)}"
    )

    # ------------------------------------------------------------------------
    # 5
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. MONOMIAL VS FALLING-BASIS COMPLEXITY")
    print("=" * 78)

    print(
        f"  F_0: monomial={complexity(F0)} "
        f"falling={complexity(fall_F0)}"
    )

    print(
        f"  F_1: monomial={complexity(F1)} "
        f"falling={complexity(fall_F1)}"
    )

    print(
        f"  H_0: monomial={complexity(H0)} "
        f"falling={complexity(H0_falling)}"
    )

    # ------------------------------------------------------------------------
    # 6
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. ROOT / FACTOR AUDIT")
    print("=" * 78)

    print(
        f"  F_0(d=0) integer roots="
        f"{integer_slice_roots(F0)}"
    )

    print(
        f"  F_1(d=0) integer roots="
        f"{integer_slice_roots(F1)}"
    )

    print(
        f"  H_0 integer roots="
        f"{boundary_source_roots(H0)}"
    )

    prop, scalar = proportional(F0, F1)

    print(
        f"  F_0/F_1 proportional={prop} "
        f"scalar={scalar}"
    )

    # ------------------------------------------------------------------------
    # 7
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. BOUNDARY-SOURCE VALUES")
    print("=" * 78)

    print(
        f"  H_0(p) = "
        f"{H0[0]} + ({H0[1]})*p"
    )

    for p in sorted(B_ODD)[:-1]:

        value = H0[0] + H0[1] * p

        print(
            f"  p={p}: H_0(p)={value}"
        )

    # ------------------------------------------------------------------------
    # 8
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. POINTWISE RECONSTRUCTION")
    print("=" * 78)

    for p in sorted(B_ODD):

        p2 = p + 2

        if p2 not in B_ODD:
            continue

        row_ok = True

        for r in range(
            max(D[p], D[p2]) + 6
        ):

            lhs = get_q(p2, r)

            rhs = evaluate_operator(
                solution,
                p,
                D[p] - r,
                r,
            )

            if lhs != rhs:
                row_ok = False

        print(
            f"  {p}->{p2}: exact={row_ok}"
        )

    # ------------------------------------------------------------------------
    # 9
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  Experiment 125 found an exact boundary-defect recurrence

      q_{p+2}(r)
        =
        F_0(p,d) q_p(r)
        +
        F_1(p,d) q_p(r+1)
        +
        H_0(p) [d=0],

  with

      d = D(p)-r.

  The present experiment converts the SAME exact operator from
  the monomial basis

      1, d, d^2, p, p*d, p^2

  to the falling basis

      1, d, d_(2), p, p*d, p_(2).

  The identities

      d^2 = d_(2) + d
      p^2 = p_(2) + p

  make this transformation exact and purely algebraic.

  Therefore any reduction in coefficient complexity is evidence
  of a genuine coordinate simplification of the already established
  operator.

  Additional tests inspect:

      * integer roots of F_0 and F_1 at d=0;
      * proportionality of F_0 and F_1;
      * roots of the boundary source H_0;
      * exact pointwise reconstruction.

  No new recurrence is fitted after the basis change.
  Everything is exact over Fraction arithmetic.
  No floating point.
  No SymPy.
  No extrapolation.
        """
    )

    # ------------------------------------------------------------------------
    # 10
    # ------------------------------------------------------------------------

    fall_roundtrip_ok = (
        fall_F0_roundtrip
        and fall_F1_roundtrip
    )

    all_ok = (
        data_exact
        and operator_exact
        and pointwise_exact
        and fall_roundtrip_ok
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(f"  data_exact={data_exact}")
    print(f"  operator_exact={operator_exact}")
    print(f"  pointwise_exact={pointwise_exact}")
    print(
        f"  falling_basis_roundtrip={fall_roundtrip_ok}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={all_ok}"
    )

    print()
    print("EXPERIMENT 126 COMPLETE")


if __name__ == "__main__":
    main()