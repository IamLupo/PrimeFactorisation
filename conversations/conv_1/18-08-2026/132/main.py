#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 132R — EXACT SYSTEM-LEVEL BASIS / DETERMINANT AUDIT
==============================================================================
"""

from fractions import Fraction
from math import gcd
from functools import reduce


# ============================================================================
# EXACT B-ODD DATA USED BY EXPERIMENTS 126-129
#
# These are the primitive integer-normalized residual rows.
#
# p=1: degree 5, six entries
# p=3: degree 4, five entries
# p=5: degree 2, three entries
# p=7: degree 0, one entry
# ============================================================================

Q = {
    1: [
        -126258696,
        -11600744760,
        2668721436,
        1764373740,
        -1338089411,
        495451247,
    ],
    3: [
        9955176,
        -1263551016,
        -152369292,
        -128667196,
        421514439,
    ],
    5: [
        -62398,
        4771718,
        16027881,
    ],
    7: [
        1,
    ],
}

D = {
    p: len(Q[p]) - 1
    for p in Q
}


# ============================================================================
# BASIC EXACT UTILITIES
# ============================================================================

def qval(p, r):
    row = Q[p]
    if 0 <= r < len(row):
        return Fraction(row[r])
    return Fraction(0)


def basis_value(name, p, d):
    if name == "one":
        return Fraction(1)

    if name == "d":
        return Fraction(d)

    if name == "d2":
        return Fraction(d * d)

    if name == "d_fall2":
        return Fraction(d * (d - 1))

    if name == "p":
        return Fraction(p)

    if name == "p2":
        return Fraction(p * p)

    if name == "p_fall2":
        return Fraction(p * (p - 1))

    if name == "p_shift":
        return Fraction(p - 1)

    if name == "p_shift_fall2":
        return Fraction((p - 1) * (p - 2))

    if name == "pd":
        return Fraction(p * d)

    if name == "p_shift_d":
        return Fraction((p - 1) * d)

    raise ValueError(f"unknown basis name: {name}")


# ============================================================================
# BASIS DICTIONARIES
# ============================================================================

BASES = {
    "monomial": {
        "F": [
            "one",
            "d",
            "d2",
            "p",
            "pd",
            "p2",
        ],
        "H": [
            "one",
            "p",
        ],
    },

    "falling": {
        "F": [
            "one",
            "d",
            "d_fall2",
            "p",
            "pd",
            "p_fall2",
        ],
        "H": [
            "one",
            "p",
        ],
    },

    "falling_reordered": {
        "F": [
            "one",
            "p",
            "p_fall2",
            "d",
            "pd",
            "d_fall2",
        ],
        "H": [
            "one",
            "p",
        ],
    },

    "shifted_falling": {
        "F": [
            "one",
            "d",
            "d_fall2",
            "p_shift",
            "p_shift_d",
            "p_shift_fall2",
        ],
        "H": [
            "one",
            "p_shift",
        ],
    },
}


# ============================================================================
# MATRIX HELPERS
# ============================================================================

def copy_matrix(A):
    return [row[:] for row in A]


def shape(A):
    return len(A), len(A[0]) if A else 0


def rank_exact(A):
    A = copy_matrix(A)

    if not A:
        return 0

    m, n = shape(A)
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

            fac = A[i][c]

            if fac == 0:
                continue

            A[i] = [
                A[i][j] - fac * A[r][j]
                for j in range(n)
            ]

        r += 1

        if r == m:
            break

    return r


def determinant(A):
    A = copy_matrix(A)
    n = len(A)

    if n == 0:
        return Fraction(1)

    if any(len(row) != n for row in A):
        raise ValueError("determinant requires a square matrix")

    det = Fraction(1)
    sign = 1

    for c in range(n):
        pivot = None

        for i in range(c, n):
            if A[i][c] != 0:
                pivot = i
                break

        if pivot is None:
            return Fraction(0)

        if pivot != c:
            A[c], A[pivot] = A[pivot], A[c]
            sign *= -1

        pv = A[c][c]
        det *= pv

        for i in range(c + 1, n):
            fac = A[i][c] / pv

            if fac == 0:
                continue

            for j in range(c + 1, n):
                A[i][j] -= fac * A[c][j]

            A[i][c] = Fraction(0)

    return det * sign


def solve_unique(A, b):
    n = len(A)

    if n == 0:
        raise ValueError("empty system")

    if any(len(row) != n for row in A):
        raise ValueError("square system required")

    M = [
        A[i][:] + [b[i]]
        for i in range(n)
    ]

    row = 0
    pivots = []

    for col in range(n):
        pivot = None

        for i in range(row, n):
            if M[i][col] != 0:
                pivot = i
                break

        if pivot is None:
            continue

        M[row], M[pivot] = M[pivot], M[row]

        pv = M[row][col]
        M[row] = [x / pv for x in M[row]]

        for i in range(n):
            if i == row:
                continue

            fac = M[i][col]

            if fac == 0:
                continue

            M[i] = [
                M[i][j] - fac * M[row][j]
                for j in range(n + 1)
            ]

        pivots.append(col)
        row += 1

    if len(pivots) != n:
        raise ValueError("system is not uniquely solvable")

    x = [Fraction(0)] * n

    for i, col in enumerate(pivots):
        x[col] = M[i][-1]

    return x


# ============================================================================
# SYSTEM CONSTRUCTION
# ============================================================================

def equation_row(p, r, basis):
    """
    Exact recurrence model:

        q_{p+2}(r)
          =
          F0(p,d) q_p(r)
          +
          F1(p,d) q_p(r+1)
          +
          H0(p) [d=0],

    d = D(p) - r.
    """

    d = D[p]

    row = []

    for name in basis["F"]:
        row.append(
            basis_value(name, p, d) * qval(p, r)
        )

    for name in basis["F"]:
        row.append(
            basis_value(name, p, d) * qval(p, r + 1)
        )

    boundary = Fraction(1) if d == 0 else Fraction(0)

    for name in basis["H"]:
        row.append(
            basis_value(name, p, d) * boundary
        )

    rhs = qval(p + 2, r)

    return row, rhs


def build_system(basis):
    """
    Exactly 14 equations:

        p=1, r=0..5 -> 6
        p=3, r=0..4 -> 5
        p=5, r=0..2 -> 3

        total = 14

    This is the same dimensionality used by the exact
    14-unknown Experiment-125/127 boundary-defect fit.
    """

    A = []
    b = []

    for p in [1, 3, 5]:
        for r in range(D[p] + 1):
            row, rhs = equation_row(p, r, basis)
            A.append(row)
            b.append(rhs)

    return A, b


# ============================================================================
# COMPLEXITY
# ============================================================================

def primitive_signature(values):
    den_lcm = 1

    for x in values:
        den_lcm = den_lcm * x.denominator // gcd(
            den_lcm,
            x.denominator,
        )

    ints = [int(x * den_lcm) for x in values]

    g = 0

    for x in ints:
        g = gcd(g, abs(x))

    if g == 0:
        return ints

    return [x // g for x in ints]


def digits(n):
    n = abs(int(n))

    if n == 0:
        return 1

    return len(str(n))


def coefficient_complexity(values):
    sig = primitive_signature(values)

    max_num = max(
        [digits(x.numerator) for x in values] or [1]
    )

    max_den = max(
        [digits(x.denominator) for x in values] or [1]
    )

    max_sig = max(
        [digits(x) for x in sig] or [1]
    )

    nonzero = sum(
        1 for x in values
        if x != 0
    )

    return (
        max_num,
        max_den,
        max_sig,
        nonzero,
    )


# ============================================================================
# CRAMER AUDIT
# ============================================================================

def cramer_solution(A, b):
    n = len(A)

    detA = determinant(A)

    if detA == 0:
        raise ValueError("singular system")

    out = []

    for col in range(n):
        M = copy_matrix(A)

        for i in range(n):
            M[i][col] = b[i]

        det_i = determinant(M)

        out.append(det_i / detA)

    return out


# ============================================================================
# BASIS ANALYSIS
# ============================================================================

def analyze(name, basis):
    A, b = build_system(basis)

    rows, cols = shape(A)
    rank = rank_exact(A)

    detA = determinant(A) if rows == cols else None

    unique = (
        rows == cols
        and rank == cols
        and detA != 0
    )

    solution = None
    reconstruction = False
    cramer_ok = False

    if unique:
        solution = solve_unique(A, b)

        reconstruction = all(
            sum(
                A[i][j] * solution[j]
                for j in range(cols)
            ) == b[i]
            for i in range(rows)
        )

        cramer = cramer_solution(A, b)

        cramer_ok = all(
            solution[i] == cramer[i]
            for i in range(cols)
        )

    common_den = 1

    if solution is not None:
        for x in solution:
            common_den = (
                common_den
                * x.denominator
                // gcd(common_den, x.denominator)
            )

    return {
        "name": name,
        "A": A,
        "b": b,
        "rows": rows,
        "cols": cols,
        "rank": rank,
        "det": detA,
        "unique": unique,
        "solution": solution,
        "reconstruction": reconstruction,
        "cramer": cramer_ok,
        "common_den": common_den,
    }


# ============================================================================
# OPERATOR FUNCTION RECOVERY
# ============================================================================

def evaluate_operator(solution, basis, p, d):
    F0 = sum(
        solution[i]
        * basis_value(basis["F"][i], p, d)
        for i in range(6)
    )

    F1 = sum(
        solution[6 + i]
        * basis_value(basis["F"][i], p, d)
        for i in range(6)
    )

    H0 = sum(
        solution[12 + i]
        * basis_value(basis["H"][i], p, d)
        for i in range(2)
    )

    return F0, F1, H0


def pointwise_operator_check(result, basis):
    if result["solution"] is None:
        return False

    sol = result["solution"]

    for p in [1, 3, 5]:

        for r in range(D[p] + 1):
            d = D[p] - r

            F0, F1, H0 = evaluate_operator(
                sol,
                basis,
                p,
                d,
            )

            lhs = qval(p + 2, r)

            rhs = (
                F0 * qval(p, r)
                + F1 * qval(p, r + 1)
                + H0 * (
                    Fraction(1)
                    if d == 0
                    else Fraction(0)
                )
            )

            if lhs != rhs:
                return False

    return True


# ============================================================================
# REPORT
# ============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 132R — EXACT SYSTEM-LEVEL BASIS / DETERMINANT AUDIT")
    print("=" * 78)
    print()

    # ------------------------------------------------------------------------
    # 1
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. EXACT DATA VALIDATION")
    print("=" * 78)

    data_exact = True

    for p in [1, 3, 5, 7]:
        degree = len(Q[p]) - 1
        expected = {
            1: 5,
            3: 4,
            5: 2,
            7: 0,
        }[p]

        ok = degree == expected

        print(
            f"  p={p}: entries={len(Q[p])} "
            f"degree={degree} expected={expected} "
            f"exact={ok}"
        )

        if not ok:
            data_exact = False

    print(f"\n  data_exact={data_exact}")
    print()

    # ------------------------------------------------------------------------
    # 2
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("2. EXACT BOUNDARY COORDINATES")
    print("=" * 78)

    for p in [1, 3, 5, 7]:
        print(
            f"  p={p}: D={D[p]} "
            f"support={list(range(D[p] + 1))}"
        )

    print()

    # ------------------------------------------------------------------------
    # 3
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("3. 14 x 14 SYSTEM VALIDATION")
    print("=" * 78)

    results = {}

    for name, basis in BASES.items():
        result = analyze(name, basis)
        results[name] = result

        print(
            f"  {name}: rows={result['rows']} "
            f"unknowns={result['cols']} "
            f"rank={result['rank']} "
            f"unique={result['unique']}"
        )

    print()

    # ------------------------------------------------------------------------
    # 4
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("4. DETERMINANT COMPARISON")
    print("=" * 78)

    for name, result in results.items():

        if result["det"] is None:
            print(f"  {name}: determinant unavailable")
            continue

        detA = result["det"]

        print(f"\n  {name}")
        print(
            f"    det_numerator_digits="
            f"{digits(detA.numerator)}"
        )
        print(
            f"    det_denominator_digits="
            f"{digits(detA.denominator)}"
        )
        print(
            f"    det_sign="
            f"{'+' if detA >= 0 else '-'}"
        )

    print()

    # ------------------------------------------------------------------------
    # 5
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("5. SOLUTION COMPLEXITY")
    print("=" * 78)

    for name, result in results.items():

        if result["solution"] is None:
            print(f"  {name}: no unique solution")
            continue

        complexity = coefficient_complexity(
            result["solution"]
        )

        print(f"\n  {name}")
        print(f"    complexity={complexity}")
        print(
            f"    common_denominator_digits="
            f"{digits(result['common_den'])}"
        )

        print(
            f"    reconstruction="
            f"{result['reconstruction']}"
        )

        print(
            f"    cramer_exact="
            f"{result['cramer']}"
        )

    print()

    # ------------------------------------------------------------------------
    # 6
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("6. POINTWISE OPERATOR RECONSTRUCTION")
    print("=" * 78)

    pointwise_ok = True

    for name, basis in BASES.items():
        result = results[name]

        ok = pointwise_operator_check(
            result,
            basis,
        )

        print(
            f"  {name}: pointwise_exact={ok}"
        )

        if not ok:
            pointwise_ok = False

    print()

    # ------------------------------------------------------------------------
    # 7
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("7. OPERATOR EQUALITY ACROSS BASES")
    print("=" * 78)

    operator_equal = True

    reference_name = None

    for name in BASES:
        if results[name]["solution"] is not None:
            reference_name = name
            break

    if reference_name is None:
        operator_equal = False
        print("  no reference operator")

    else:
        reference_basis = BASES[reference_name]
        reference_solution = results[
            reference_name
        ]["solution"]

        for name, basis in BASES.items():

            if results[name]["solution"] is None:
                operator_equal = False
                print(
                    f"  {reference_name} vs {name}: "
                    f"unavailable"
                )
                continue

            other_solution = results[name]["solution"]

            equivalent = True

            for p in [1, 3, 5]:

                for r in range(D[p] + 1):
                    d = D[p] - r

                    F0a, F1a, H0a = evaluate_operator(
                        reference_solution,
                        reference_basis,
                        p,
                        d,
                    )

                    F0b, F1b, H0b = evaluate_operator(
                        other_solution,
                        basis,
                        p,
                        d,
                    )

                    if (
                        F0a != F0b
                        or F1a != F1b
                        or H0a != H0b
                    ):
                        equivalent = False
                        break

                if not equivalent:
                    break

            print(
                f"  {reference_name} vs {name}: "
                f"same_operator={equivalent}"
            )

            if not equivalent:
                operator_equal = False

    print()

    # ------------------------------------------------------------------------
    # 8
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("8. FALLING-BASIS IDENTITIES")
    print("=" * 78)

    identity_ok = True

    for x in range(0, 10):

        if x * x != x * (x - 1) + x:
            identity_ok = False

    print(
        "  d^2 = d_(2) + d : "
        f"{identity_ok}"
    )

    print(
        "  p^2 = p_(2) + p : "
        f"{identity_ok}"
    )

    print()

    # ------------------------------------------------------------------------
    # 9
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  The previous Experiment 132 mixed two different data conventions.

  This corrected audit uses exactly the primitive integer-normalized
  B-odd rows from the hardened Experiment-129 data:

      p=1 : degree 5
      p=3 : degree 4
      p=5 : degree 2
      p=7 : degree 0

  Hence the boundary coordinates are

      D(1)=5,
      D(3)=4,
      D(5)=2,
      D(7)=0.

  The 14 equations are exactly

      6 + 5 + 3 = 14,

  matching the 14 unknowns in the boundary-defect operator:

      F_0 : 6 coefficients,
      F_1 : 6 coefficients,
      H_0 : 2 coefficients.

  The experiment rebuilds the linear system directly in several
  coordinate bases before solving it.

  Therefore this is a genuine system-level test of the determinant
  inflation hypothesis.

  The relevant question is not whether the final operator can be
  rewritten more elegantly, but whether the interpolation matrix
  itself becomes algebraically simpler before solving.

  In particular:

      lower determinant size
      lower denominator size
      lower primitive-signature size
      fewer nonzero coordinates

  would support the hypothesis that the enormous rational coefficients
  are partly caused by a poor coordinate basis.

  No connection to the original (p,q)-kernel is assumed here.

  Everything is exact Fraction arithmetic.
  No floating point.
  No SymPy.
  No extrapolation.
        """.strip()
    )

    print()

    # ------------------------------------------------------------------------
    # 10
    # ------------------------------------------------------------------------

    all_unique = all(
        results[name]["unique"]
        for name in BASES
    )

    all_reconstruction = all(
        results[name]["reconstruction"]
        for name in BASES
    )

    all_cramer = all(
        results[name]["cramer"]
        for name in BASES
    )

    failures = 0

    if not data_exact:
        failures += 1

    if not all_unique:
        failures += 1

    if not all_reconstruction:
        failures += 1

    if not all_cramer:
        failures += 1

    if not pointwise_ok:
        failures += 1

    if not operator_equal:
        failures += 1

    if not identity_ok:
        failures += 1

    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  data_exact={data_exact}"
    )
    print(
        f"  all_basis_systems_unique={all_unique}"
    )
    print(
        f"  all_reconstructions_exact={all_reconstruction}"
    )
    print(
        f"  all_cramer_checks_exact={all_cramer}"
    )
    print(
        f"  pointwise_operator_exact={pointwise_ok}"
    )
    print(
        f"  cross_basis_operator_equivalence={operator_equal}"
    )
    print(
        f"  falling_basis_identity={identity_ok}"
    )
    print(
        f"  failures={failures}"
    )
    print(
        f"  ALL BASIC CHECKS PASS={failures == 0}"
    )

    print()
    print("EXPERIMENT 132R COMPLETE")


if __name__ == "__main__":
    main()