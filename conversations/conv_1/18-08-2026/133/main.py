#!/usr/bin/env python3

from fractions import Fraction
from math import gcd


# ============================================================================
# EXPERIMENT 132S
# EXACT SYSTEM-LEVEL FALLING-BASIS / DETERMINANT AUDIT
#
# CRITICAL:
# Use the actual rational q[p,r] rows.
# Do NOT replace them by primitive integer signatures.
# ============================================================================


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
    p: len(Q[p]) - 1
    for p in Q
}


# ============================================================================
# BASIC EXACT UTILITIES
# ============================================================================

def qval(p, r):
    if 0 <= r < len(Q[p]):
        return Q[p][r]
    return Fraction(0)


def degree_row(p):
    return len(Q[p]) - 1


def falling2(x):
    return x * (x - 1)


def basis_value(name, p, d):
    if name == "one":
        return Fraction(1)
    if name == "d":
        return Fraction(d)
    if name == "d2":
        return Fraction(d * d)
    if name == "p":
        return Fraction(p)
    if name == "pd":
        return Fraction(p * d)
    if name == "p2":
        return Fraction(p * p)

    if name == "d_fall2":
        return Fraction(falling2(d))
    if name == "p_fall2":
        return Fraction(falling2(p))
    if name == "p_shift":
        return Fraction(p - 1)
    if name == "p_shift_d":
        return Fraction((p - 1) * d)
    if name == "p_shift_fall2":
        return Fraction(falling2(p - 1))

    raise ValueError("Unknown basis name: " + str(name))


# ============================================================================
# RATIONAL MATRIX ROUTINES
# ============================================================================

def matrix_copy(A):
    return [row[:] for row in A]


def rank_exact(A):
    if not A:
        return 0

    M = matrix_copy(A)
    m = len(M)
    n = len(M[0])
    rank = 0

    for col in range(n):
        pivot = None

        for row in range(rank, m):
            if M[row][col] != 0:
                pivot = row
                break

        if pivot is None:
            continue

        if pivot != rank:
            M[rank], M[pivot] = M[pivot], M[rank]

        pivot_value = M[rank][col]
        M[rank] = [
            x / pivot_value
            for x in M[rank]
        ]

        for row in range(m):
            if row == rank:
                continue

            factor = M[row][col]

            if factor == 0:
                continue

            M[row] = [
                M[row][j] - factor * M[rank][j]
                for j in range(n)
            ]

        rank += 1

        if rank == m:
            break

    return rank


def solve_unique(A, b):
    n = len(A)

    if n == 0:
        raise ValueError("Empty system.")

    if any(len(row) != n for row in A):
        raise ValueError("System is not square.")

    M = [
        A[i][:] + [b[i]]
        for i in range(n)
    ]

    pivot_row = 0
    pivots = []

    for col in range(n):

        pivot = None

        for row in range(pivot_row, n):
            if M[row][col] != 0:
                pivot = row
                break

        if pivot is None:
            continue

        if pivot != pivot_row:
            M[pivot_row], M[pivot] = (
                M[pivot],
                M[pivot_row],
            )

        pv = M[pivot_row][col]

        M[pivot_row] = [
            x / pv
            for x in M[pivot_row]
        ]

        for row in range(n):

            if row == pivot_row:
                continue

            factor = M[row][col]

            if factor == 0:
                continue

            M[row] = [
                M[row][j]
                - factor * M[pivot_row][j]
                for j in range(n + 1)
            ]

        pivots.append(col)
        pivot_row += 1

    if len(pivots) != n:
        raise ValueError("System is not uniquely solvable.")

    result = [Fraction(0)] * n

    for row, col in enumerate(pivots):
        result[col] = M[row][-1]

    return result


def determinant(A):
    n = len(A)

    if n == 0:
        return Fraction(1)

    if any(len(row) != n for row in A):
        raise ValueError("Determinant requires a square matrix.")

    M = matrix_copy(A)
    det = Fraction(1)
    sign = 1

    for col in range(n):

        pivot = None

        for row in range(col, n):
            if M[row][col] != 0:
                pivot = row
                break

        if pivot is None:
            return Fraction(0)

        if pivot != col:
            M[col], M[pivot] = (
                M[pivot],
                M[col],
            )
            sign *= -1

        pv = M[col][col]
        det *= pv

        for row in range(col + 1, n):

            factor = M[row][col] / pv

            if factor == 0:
                continue

            for j in range(col + 1, n):
                M[row][j] -= factor * M[col][j]

            M[row][col] = Fraction(0)

    return det * sign


# ============================================================================
# EXACT BASIS DEFINITIONS
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
# SYSTEM
# ============================================================================

def build_system(basis):

    A = []
    b = []

    # 6 + 5 + 3 = 14 equations.
    #
    # p=1 -> r=0,...,5
    # p=3 -> r=0,...,4
    # p=5 -> r=0,...,2
    #
    # The p=7 row is the terminal target.

    for p in [1, 3, 5]:

        for r in range(D[p] + 1):

            d = D[p] - r

            row = []

            # F_0(p,d) q_p(r)
            for name in basis["F"]:
                row.append(
                    basis_value(name, p, d)
                    * qval(p, r)
                )

            # F_1(p,d) q_p(r+1)
            for name in basis["F"]:
                row.append(
                    basis_value(name, p, d)
                    * qval(p, r + 1)
                )

            # H_0(p) [d=0]
            indicator = (
                Fraction(1)
                if d == 0
                else Fraction(0)
            )

            for name in basis["H"]:
                row.append(
                    basis_value(name, p, d)
                    * indicator
                )

            A.append(row)
            b.append(qval(p + 2, r))

    return A, b


# ============================================================================
# INTEGER SIGNATURE
# ============================================================================

def primitive_signature(values):

    lcm_den = 1

    for x in values:
        a = x.denominator
        g = gcd(lcm_den, a)
        lcm_den = lcm_den // g * a

    ints = [
        int(x * lcm_den)
        for x in values
    ]

    common = 0

    for x in ints:
        common = gcd(common, abs(x))

    if common == 0:
        return ints

    return [
        x // common
        for x in ints
    ]


def digit_count(x):

    x = abs(int(x))

    if x == 0:
        return 1

    return len(str(x))


def complexity(values):

    sig = primitive_signature(values)

    max_num = max(
        [digit_count(x.numerator) for x in values]
        or [1]
    )

    max_den = max(
        [digit_count(x.denominator) for x in values]
        or [1]
    )

    max_sig = max(
        [digit_count(x) for x in sig]
        or [1]
    )

    nonzero = sum(
        1
        for x in values
        if x != 0
    )

    return (
        max_num,
        max_den,
        max_sig,
        nonzero,
    )


# ============================================================================
# OPERATOR EVALUATION
# ============================================================================

def evaluate_operator(solution, basis, p, d):

    F0 = Fraction(0)
    F1 = Fraction(0)
    H0 = Fraction(0)

    for i, name in enumerate(basis["F"]):

        v = basis_value(
            name,
            p,
            d,
        )

        F0 += solution[i] * v
        F1 += solution[6 + i] * v

    for i, name in enumerate(basis["H"]):

        H0 += (
            solution[12 + i]
            * basis_value(name, p, d)
        )

    return F0, F1, H0


def pointwise_check(solution, basis):

    if solution is None:
        return False

    for p in [1, 3, 5]:

        for r in range(D[p] + 1):

            d = D[p] - r

            F0, F1, H0 = evaluate_operator(
                solution,
                basis,
                p,
                d,
            )

            rhs = (
                F0 * qval(p, r)
                +
                F1 * qval(p, r + 1)
                +
                H0
                * (
                    Fraction(1)
                    if d == 0
                    else Fraction(0)
                )
            )

            lhs = qval(p + 2, r)

            if lhs != rhs:
                return False

    return True


# ============================================================================
# CRAMER CHECK
# ============================================================================

def cramer_check(A, b, solution):

    if solution is None:
        return False

    detA = determinant(A)

    if detA == 0:
        return False

    for col in range(len(A)):

        M = matrix_copy(A)

        for row in range(len(A)):
            M[row][col] = b[row]

        if determinant(M) != solution[col] * detA:
            return False

    return True


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 132S — EXACT SYSTEM-LEVEL")
    print("FALLING-BASIS / DETERMINANT AUDIT")
    print("=" * 78)
    print()

    # ------------------------------------------------------------------------
    # 1
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. EXACT DATA VALIDATION")
    print("=" * 78)

    expected_degree = {
        1: 5,
        3: 4,
        5: 2,
        7: 0,
    }

    data_exact = True

    for p in [1, 3, 5, 7]:

        deg = degree_row(p)

        ok = deg == expected_degree[p]

        print(
            f"  p={p}: entries={len(Q[p])} "
            f"degree={deg} "
            f"expected={expected_degree[p]} "
            f"exact={ok}"
        )

        if not ok:
            data_exact = False

    print()
    print(
        f"  data_exact={data_exact}"
    )
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
    print("3. EXACT 14x14 SYSTEM RANK")
    print("=" * 78)

    results = {}

    for name, basis in BASES.items():

        A, b = build_system(basis)

        rankA = rank_exact(A)
        augmented = [
            A[i][:] + [b[i]]
            for i in range(len(A))
        ]

        rankAug = rank_exact(augmented)

        detA = determinant(A)

        unique = (
            len(A) == 14
            and len(A[0]) == 14
            and rankA == 14
            and rankAug == 14
            and detA != 0
        )

        solution = None

        if unique:
            solution = solve_unique(A, b)

        results[name] = {
            "A": A,
            "b": b,
            "rank": rankA,
            "augmented_rank": rankAug,
            "det": detA,
            "unique": unique,
            "solution": solution,
        }

        print(
            f"  {name}: "
            f"rows={len(A)} "
            f"unknowns={len(A[0])} "
            f"rank={rankA} "
            f"augmented_rank={rankAug} "
            f"unique={unique}"
        )

    print()

    # ------------------------------------------------------------------------
    # 4
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("4. DETERMINANT SIZE")
    print("=" * 78)

    for name, result in results.items():

        detA = result["det"]

        print(f"\n  {name}")
        print(
            f"    numerator_digits="
            f"{digit_count(detA.numerator)}"
        )
        print(
            f"    denominator_digits="
            f"{digit_count(detA.denominator)}"
        )
        print(
            f"    sign="
            f"{'+' if detA > 0 else '-'}"
            if detA != 0
            else
            "    sign=0"
        )

    print()

    # ------------------------------------------------------------------------
    # 5
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("5. SOLUTION COMPLEXITY")
    print("=" * 78)

    for name, result in results.items():

        sol = result["solution"]

        if sol is None:
            print(
                f"  {name}: no unique exact solution"
            )
            continue

        sig = primitive_signature(sol)

        print(f"\n  {name}")
        print(
            f"    complexity={complexity(sol)}"
        )
        print(
            f"    max_num_digits="
            f"{max(digit_count(x.numerator) for x in sol)}"
        )
        print(
            f"    max_den_digits="
            f"{max(digit_count(x.denominator) for x in sol)}"
        )
        print(
            f"    primitive_signature={sig}"
        )

    print()

    # ------------------------------------------------------------------------
    # 6
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("6. POINTWISE OPERATOR RECONSTRUCTION")
    print("=" * 78)

    pointwise_all = True

    for name, basis in BASES.items():

        sol = results[name]["solution"]

        ok = pointwise_check(
            sol,
            basis,
        )

        print(
            f"  {name}: pointwise_exact={ok}"
        )

        if not ok:
            pointwise_all = False

    print()

    # ------------------------------------------------------------------------
    # 7
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("7. EXACT CRAMER AUDIT")
    print("=" * 78)

    cramer_all = True

    for name, result in results.items():

        ok = cramer_check(
            result["A"],
            result["b"],
            result["solution"],
        )

        print(
            f"  {name}: "
            f"cramer_exact={ok}"
        )

        if not ok:
            cramer_all = False

    print()

    # ------------------------------------------------------------------------
    # 8
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("8. FALLING-BASIS IDENTITIES")
    print("=" * 78)

    d_ok = all(
        d * d == d * (d - 1) + d
        for d in range(0, 20)
    )

    p_ok = all(
        p * p == p * (p - 1) + p
        for p in range(0, 20)
    )

    print(
        f"  d^2 = d_(2) + d : {d_ok}"
    )
    print(
        f"  p^2 = p_(2) + p : {p_ok}"
    )

    print()

    # ------------------------------------------------------------------------
    # 9
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("9. CRITICAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  IMPORTANT DATA DISTINCTION

  The primitive integer signatures printed by Experiments 129-131
  are NOT the same objects as the actual q[p,r] coefficients.

  Independent row rescaling is harmless for some rank questions,
  but it is NOT harmless for an operator equation of the form

      q_{p+2}(r)
        =
        F_0(p,d) q_p(r)
        +
        F_1(p,d) q_p(r+1)
        +
        H_0(p) [d=0].

  Therefore the recurrence must be reconstructed from the exact
  rational q[p,r] rows, not from their primitive signatures.

  This script uses the exact rational rows from Experiment 115:

      q_1 : 6 entries
      q_3 : 5 entries
      q_5 : 3 entries
      q_7 : 1 entry.

  Consequently the exact system again has

      6 + 5 + 3 = 14 equations

  and

      6 + 6 + 2 = 14 unknowns.

  Only after this correction does it make sense to investigate
  determinant inflation, Cramer's-rule denominator growth, and
  basis-dependent coefficient complexity.

  No connection to the original (p,q)-kernel is asserted here.

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
        result["unique"]
        for result in results.values()
    )

    all_pointwise = pointwise_all

    failures = 0

    if not data_exact:
        failures += 1

    if not all_unique:
        failures += 1

    if not all_pointwise:
        failures += 1

    if not cramer_all:
        failures += 1

    if not d_ok or not p_ok:
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
        f"  pointwise_operator_exact={all_pointwise}"
    )
    print(
        f"  cramer_audit_exact={cramer_all}"
    )
    print(
        f"  falling_basis_identity="
        f"{d_ok and p_ok}"
    )
    print(
        f"  failures={failures}"
    )
    print(
        f"  ALL BASIC CHECKS PASS={failures == 0}"
    )

    print()
    print("EXPERIMENT 132S COMPLETE")


if __name__ == "__main__":
    main()

