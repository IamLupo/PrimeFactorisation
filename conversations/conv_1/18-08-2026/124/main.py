from fractions import Fraction
from math import gcd
from itertools import product


# ============================================================================
# EXPERIMENT 124 — EXACT BOUNDARY-AWARE p/r SHIFT-OPERATOR AUDIT
# ============================================================================
#
# Goal:
#   Experiment 123 found no UNIQUE exact operator of total degree <= 2
#   in the raw variables (p,r).
#
#   The next natural coordinate is the finite-support boundary distance
#
#       d = D(p) - r,
#
#   where D(p) is the observed last nonzero falling-basis index.
#
#   We search
#
#       q_{p+2}(r)
#         = sum_{a=0}^m F_a(p,d) q_p(r+a)
#
#   with
#
#       F_a(p,d)
#
#   polynomial in (p,d) of bounded total degree.
#
#   This explicitly encodes the finite-support boundary instead of
#   treating r as an unbounded coordinate.
#
#   Arithmetic:
#       exact Fraction only
#       no floating point
#       no SymPy
#       no extrapolation
#
# ============================================================================


# ----------------------------------------------------------------------------
# 1. EXACT DATA
# ----------------------------------------------------------------------------

# B-odd corrected second-layer rows from Experiment 115.
#
# q_p[r] means the coefficient of r_(r) in Q_p(r).
#
# p = 1,3,5,7
#
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

# Last nonzero r-index D(p).
D = {p: len(row) - 1 for p, row in B_ODD.items()}


# ----------------------------------------------------------------------------
# 2. OPTIONAL OTHER SECTORS
# ----------------------------------------------------------------------------
#
# The experiment first focuses on the exceptional B-odd sector.
# These exact rows are included so the same operator can later be tested
# against the other sectors without rewriting the infrastructure.
# ----------------------------------------------------------------------------

A_EVEN = {
    0: [
        Fraction(-12879, 1),
        Fraction(55745625600, 46735315200),
        Fraction(-30144441600, 46735315200),
        Fraction(2714342400, 46735315200),
        Fraction(1570086000, 46735315200),
        Fraction(-613274346, 46735315200),
        Fraction(102402481, 46735315200),
    ],
    2: [
        Fraction(2797337, 11520),
        Fraction(44270519328, 46735315200),
        Fraction(-19166962968, 46735315200),
        Fraction(-1132410616, 46735315200),
        Fraction(2461903867, 46735315200),
        Fraction(-783039371, 46735315200),
    ],
    4: [
        Fraction(-2083937, 921600),
        Fraction(1980796104, 46735315200),
        Fraction(-1311194478, 46735315200),
        Fraction(-477293762, 46735315200),
        Fraction(727023095, 46735315200),
    ],
    6: [
        Fraction(85591, 22118400),
        Fraction(6972364, 46735315200),
        Fraction(-24042833, 46735315200),
    ],
    8: [
        Fraction(-4913, 353894400),
    ],
}

# Keep the main experiment B-odd.  The extra data are retained only for
# optional cross-sector boundary testing.


# ----------------------------------------------------------------------------
# 3. BASIC EXACT ARITHMETIC / LINEAR ALGEBRA
# ----------------------------------------------------------------------------

def frac_str(x):
    return str(Fraction(x))


def matrix_rank(A):
    """Exact rank over Q."""
    if not A:
        return 0

    M = [list(map(Fraction, row)) for row in A]
    rows = len(M)
    cols = len(M[0]) if rows else 0

    rank = 0
    col = 0

    while rank < rows and col < cols:
        pivot = None
        for r in range(rank, rows):
            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            col += 1
            continue

        if pivot != rank:
            M[rank], M[pivot] = M[pivot], M[rank]

        pv = M[rank][col]
        for j in range(col, cols):
            M[rank][j] /= pv

        for r in range(rows):
            if r == rank:
                continue
            if M[r][col] == 0:
                continue

            factor = M[r][col]
            for j in range(col, cols):
                M[r][j] -= factor * M[rank][j]

        rank += 1
        col += 1

    return rank


def rref(A):
    """Exact RREF over Q."""
    if not A:
        return [], []

    M = [list(map(Fraction, row)) for row in A]
    rows = len(M)
    cols = len(M[0]) if rows else 0

    pivots = []
    rank = 0
    col = 0

    while rank < rows and col < cols:
        pivot = None

        for r in range(rank, rows):
            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            col += 1
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
        col += 1

    return M, pivots


def solve_exact(A, b):
    """
    Exact linear system.
    Returns:
        None                     inconsistent
        []                       no unknowns
        list of one solution     unique
    """
    if not A:
        return []

    if len(A) != len(b):
        raise ValueError("row mismatch")

    aug = [
        list(map(Fraction, A[i])) + [Fraction(b[i])]
        for i in range(len(A))
    ]

    R, pivots = rref(aug)

    n = len(A[0])

    # inconsistency
    for row in R:
        if all(row[j] == 0 for j in range(n)) and row[n] != 0:
            return None

    if len(pivots) < n:
        return None

    sol = [Fraction(0) for _ in range(n)]

    for i, p in enumerate(pivots):
        sol[p] = R[i][n]

    return sol


# ----------------------------------------------------------------------------
# 4. POLYNOMIAL BASIS IN (p,d)
# ----------------------------------------------------------------------------

def basis_exponents(total_degree):
    out = []
    for i in range(total_degree + 1):
        for j in range(total_degree + 1 - i):
            out.append((i, j))
    return out


def eval_basis(p, d, degree):
    """
    Monomial basis:
        p^i d^j,  i+j <= degree
    """
    exps = basis_exponents(degree)
    return [Fraction(p ** i) * Fraction(d ** j) for i, j in exps]


def basis_to_string(coeffs, degree, name="F"):
    exps = basis_exponents(degree)
    terms = []

    for c, (i, j) in zip(coeffs, exps):
        if c == 0:
            continue

        if i == 0 and j == 0:
            mon = ""
        elif i == 0:
            mon = "d" if j == 1 else f"d^{j}"
        elif j == 0:
            mon = "p" if i == 1 else f"p^{i}"
        else:
            ppart = "p" if i == 1 else f"p^{i}"
            dpart = "d" if j == 1 else f"d^{j}"
            mon = f"{ppart}*{dpart}"

        if mon:
            terms.append(f"({c})*{mon}")
        else:
            terms.append(f"({c})")

    if not terms:
        return "0"

    return " + ".join(terms)


# ----------------------------------------------------------------------------
# 5. GENERAL BOUNDARY-AWARE OPERATOR SYSTEM
# ----------------------------------------------------------------------------

def get_q(row, r):
    if 0 <= r < len(row):
        return row[r]
    return Fraction(0)


def build_operator_system(rows, Dmap, order, degree, step=2):
    """
    Search

        q[p+step, r]
          = sum_a F_a(p,d) q[p, r+a]

    where d = D(p)-r,

    and every F_a is a polynomial in p,d of total degree <= degree.

    Unknown ordering:
        (a, basis_index)

    We only generate equations where the target p+step exists.

    Boundary-zero equations are included automatically because get_q()
    returns 0 outside the finite support.
    """
    ps = sorted(rows)
    transitions = [(p, p + step) for p in ps if p + step in rows]

    exps = basis_exponents(degree)
    nb = len(exps)
    n_unknowns = order + 1
    nvars = n_unknowns * nb

    A = []
    b = []

    for p, p2 in transitions:
        Dp = Dmap[p]

        # Include all integer r values that can influence either side.
        max_r = max(
            Dp,
            Dmap[p2],
            Dp + order,
            Dmap[p2] + order,
        )

        for r in range(0, max_r + 1):
            d = Dp - r

            basis = eval_basis(p, d, degree)

            row_eq = [Fraction(0) for _ in range(nvars)]

            for a in range(order + 1):
                qval = get_q(rows[p], r + a)

                if qval == 0:
                    continue

                offset = a * nb

                for j, value in enumerate(basis):
                    row_eq[offset + j] += qval * value

            target = get_q(rows[p2], r)

            # Equation:
            # sum_a F_a * q_p(r+a) - q_p2(r) = 0
            A.append(row_eq)
            b.append(target)

    return A, b, nb, exps


def operator_search(rows, Dmap, order, degree, step=2):
    A, b, nb, exps = build_operator_system(
        rows, Dmap, order, degree, step
    )

    if not A:
        return {
            "consistent": False,
            "unique": False,
            "rank": 0,
            "unknowns": order * nb,
            "nullity": None,
            "solution": None,
            "basis": exps,
        }

    # We solve:
    # A*x = b
    sol = solve_exact(A, b)
    rankA = matrix_rank(A)
    rankAug = matrix_rank([
        row + [b[i]]
        for i, row in enumerate(A)
    ])

    nvars = order * nb

    if sol is None:
        return {
            "consistent": False,
            "unique": False,
            "rank": rankA,
            "unknowns": nvars,
            "nullity": None,
            "solution": None,
            "basis": exps,
        }

    nullity = nvars - rankA
    unique = nullity == 0

    return {
        "consistent": rankA == rankAug,
        "unique": unique,
        "rank": rankA,
        "unknowns": nvars,
        "nullity": nullity,
        "solution": sol,
        "basis": exps,
    }


# ----------------------------------------------------------------------------
# 6. SPECIAL BOUNDARY CONSTRAINTS
# ----------------------------------------------------------------------------

def coefficient_block(solution, a, degree):
    nb = len(basis_exponents(degree))
    start = a * nb
    return solution[start:start + nb]


def evaluate_block(block, p, d, degree):
    basis = eval_basis(p, d, degree)
    return sum(c * v for c, v in zip(block, basis))


def boundary_zero_residuals(rows, Dmap, solution, order, degree, step=2):
    """
    Check equations at boundary points separately.

    For a finite-support row, a shift a that goes beyond the previous
    support must be multiplied by a coefficient that vanishes whenever
    that shifted term is impossible.

    We report the exact values of F_a at those boundary configurations.
    """
    checks = []

    blocks = [
        coefficient_block(solution, a, degree)
        for a in range(order + 1)
    ]

    for p in sorted(rows):
        p2 = p + step

        if p2 not in rows:
            continue

        Dp = Dmap[p]
        Dp2 = Dmap[p2]

        # At each r >= 0, terms with r+a > Dp are structurally absent.
        for r in range(0, Dp2 + 1):
            d = Dp - r

            for a in range(order + 1):
                source_exists = (r + a) <= Dp
                if source_exists:
                    continue

                val = evaluate_block(
                    blocks[a],
                    p,
                    d,
                    degree,
                )

                checks.append(
                    (p, p2, r, a, val)
                )

    return checks


def all_boundary_constraints_zero(checks):
    return all(item[-1] == 0 for item in checks)


# ----------------------------------------------------------------------------
# 7. DIVISIBILITY / FACTOR TESTS IN d
# ----------------------------------------------------------------------------

def univariate_integer_roots_in_range(values, lo, hi):
    """
    Evaluate a polynomial represented by dense ascending coefficients
    at integer points. Used only for small exact factor checks.
    """
    hits = []
    for x in range(lo, hi + 1):
        total = Fraction(0)
        for i, c in enumerate(values):
            total += c * Fraction(x ** i)
        if total == 0:
            hits.append(x)
    return hits


def block_values_at_p(block, degree, p, d_values):
    return [
        evaluate_block(block, p, d, degree)
        for d in d_values
    ]


# ----------------------------------------------------------------------------
# 8. OPERATOR POINTWISE VERIFICATION
# ----------------------------------------------------------------------------

def verify_operator(rows, Dmap, solution, order, degree, step=2):
    blocks = [
        coefficient_block(solution, a, degree)
        for a in range(order + 1)
    ]

    failures = []

    ps = sorted(rows)
    for p in ps:
        p2 = p + step
        if p2 not in rows:
            continue

        Dp2 = Dmap[p2]

        for r in range(0, Dp2 + 1):
            d = Dmap[p] - r

            lhs = get_q(rows[p2], r)

            rhs = Fraction(0)
            for a in range(order + 1):
                rhs += (
                    evaluate_block(
                        blocks[a],
                        p,
                        d,
                        degree,
                    )
                    * get_q(rows[p], r + a)
                )

            if lhs != rhs:
                failures.append(
                    (p, p2, r, lhs, rhs)
                )

    return failures


# ----------------------------------------------------------------------------
# 9. SMALL FACTORIZATION TEST
# ----------------------------------------------------------------------------

def block_integer_signature(block):
    """
    Convert rational coefficient block to a primitive integer signature.
    """
    if not block:
        return []

    den_lcm = 1

    for c in block:
        den = c.denominator
        den_lcm = den_lcm * den // gcd(den_lcm, den)

    ints = [int(c * den_lcm) for c in block]

    g = 0
    for n in ints:
        g = gcd(g, abs(n))

    if g == 0:
        return [0 for _ in ints]

    ints = [n // g for n in ints]

    for n in ints:
        if n != 0:
            if n < 0:
                ints = [-x for x in ints]
            break

    return ints


def print_solution_blocks(solution, order, degree, title="operator"):
    print(f"  {title}:")
    for a in range(order + 1):
        block = coefficient_block(solution, a, degree)
        print(
            f"    F_{a}(p,d) = {basis_to_string(block, degree)}"
        )
        print(
            f"      primitive_signature = "
            f"{block_integer_signature(block)}"
        )


# ----------------------------------------------------------------------------
# 10. SEARCH OVER THE BOUNDARY-AWARE MODELS
# ----------------------------------------------------------------------------

def systematic_search(rows, Dmap):
    print("=" * 78)
    print("BOUNDARY-AWARE OPERATOR SEARCH")
    print("=" * 78)

    records = []

    for order in (1, 2, 3):
        for degree in (0, 1, 2, 3):

            result = operator_search(
                rows,
                Dmap,
                order,
                degree,
                step=2,
            )

            records.append(
                (order, degree, result)
            )

            print(
                f"  order={order} degree={degree}: "
                f"consistent={result['consistent']} "
                f"rank={result['rank']} "
                f"unknowns={result['unknowns']} "
                f"nullity={result['nullity']}"
            )

            if result["consistent"] and result["solution"] is not None:
                failures = verify_operator(
                    rows,
                    Dmap,
                    result["solution"],
                    order,
                    degree,
                )

                boundary_checks = boundary_zero_residuals(
                    rows,
                    Dmap,
                    result["solution"],
                    order,
                    degree,
                )

                boundary_ok = all_boundary_constraints_zero(
                    boundary_checks
                )

                print(
                    f"    pointwise_exact={len(failures) == 0}"
                )
                print(
                    f"    boundary_zero_constraints="
                    f"{len(boundary_checks)}"
                )
                print(
                    f"    boundary_zero_exact={boundary_ok}"
                )

                if result["unique"]:
                    print("    UNIQUE-EXACT")

    return records


def select_best_unique(records):
    candidates = []

    for order, degree, result in records:
        if not (
            result["consistent"]
            and result["unique"]
            and result["solution"] is not None
        ):
            continue

        candidates.append(
            (order, degree, result)
        )

    if not candidates:
        return None

    # Prefer:
    #   1. smaller order
    #   2. smaller degree
    return min(
        candidates,
        key=lambda item: (item[0], item[1])
    )


# ----------------------------------------------------------------------------
# 11. BOUNDARY FACTOR AUDIT
# ----------------------------------------------------------------------------

def boundary_factor_audit(rows, Dmap, solution, order, degree):
    print()
    print("=" * 78)
    print("BOUNDARY-FACTOR AUDIT")
    print("=" * 78)

    blocks = [
        coefficient_block(solution, a, degree)
        for a in range(order + 1)
    ]

    any_factor = False

    for a, block in enumerate(blocks):

        print(f"  F_{a}(p,d)")

        for p in sorted(rows):
            max_d = Dmap[p]

            vals = [
                evaluate_block(
                    block,
                    p,
                    d,
                    degree,
                )
                for d in range(0, max_d + 1)
            ]

            zero_positions = [
                d
                for d, value in enumerate(vals)
                if value == 0
            ]

            if zero_positions:
                any_factor = True

            print(
                f"    p={p}: "
                f"d_zeros={zero_positions}"
            )

    print(
        f"  any_exact_boundary_zero={any_factor}"
    )

    return any_factor


# ----------------------------------------------------------------------------
# 12. CROSS-SECTOR TEST
# ----------------------------------------------------------------------------

def cross_sector_test():

    print()
    print("=" * 78)
    print("CROSS-SECTOR BOUNDARY-AWARE TRANSFER TEST")
    print("=" * 78)

    # Only sectors whose exact second-layer rows are available here.
    #
    # A-even is deliberately included as a contrasting sector.
    sectors = {
        "A-even": (
            A_EVEN,
            {p: len(row) - 1 for p, row in A_EVEN.items()},
        ),
        "B-odd": (
            B_ODD,
            D,
        ),
    }

    for name, (rows, dmap) in sectors.items():

        print(f"  {name}")

        found = []

        for order in (1, 2, 3):
            for degree in (0, 1, 2, 3):

                res = operator_search(
                    rows,
                    dmap,
                    order,
                    degree,
                    step=2,
                )

                if (
                    res["consistent"]
                    and res["unique"]
                    and res["solution"] is not None
                ):
                    fails = verify_operator(
                        rows,
                        dmap,
                        res["solution"],
                        order,
                        degree,
                    )

                    if not fails:
                        found.append(
                            (order, degree)
                        )

        print(
            f"    unique_exact_boundary_models={found}"
        )


# ----------------------------------------------------------------------------
# 13. DATA SANITY
# ----------------------------------------------------------------------------

def sanity_check(rows, Dmap):

    print()
    print("=" * 78)
    print("1. EXACT DATA VALIDATION")
    print("=" * 78)

    ok = True

    for p in sorted(rows):

        row = rows[p]

        expected_degree = Dmap[p]
        actual_degree = len(row) - 1

        row_ok = actual_degree == expected_degree

        print(
            f"  p={p}: "
            f"entries={len(row)} "
            f"degree={actual_degree} "
            f"expected={expected_degree} "
            f"exact={row_ok}"
        )

        ok = ok and row_ok

    for p in sorted(rows):
        p2 = p + 2

        if p2 not in rows:
            continue

        if len(rows[p2]) > len(rows[p]) + 1:
            # The next row cannot have degree growth beyond the
            # available source neighborhood for a small forward operator.
            print(
                f"  transition {p}->{p2}: "
                f"support_growth_warning=True"
            )

    print(f"  data_exact={ok}")
    return ok


# ----------------------------------------------------------------------------
# 14. MAIN
# ----------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 124 — EXACT BOUNDARY-AWARE p/r SHIFT-OPERATOR AUDIT")
    print("=" * 78)
    print()

    data_ok = sanity_check(
        B_ODD,
        D,
    )

    print()
    print("=" * 78)
    print("2. OBSERVED BOUNDARY COORDINATES")
    print("=" * 78)

    for p in sorted(B_ODD):
        print(
            f"  p={p}: "
            f"D={D[p]} "
            f"support={list(range(D[p] + 1))}"
        )

    print()
    print("=" * 78)
    print("3. SYSTEMATIC BOUNDARY-AWARE SEARCH")
    print("=" * 78)

    records = systematic_search(
        B_ODD,
        D,
    )

    best = select_best_unique(records)

    print()
    print("=" * 78)
    print("4. LOWEST-COMPLEXITY UNIQUE EXACT MODEL")
    print("=" * 78)

    best_ok = False

    if best is None:
        print("  NONE")
    else:
        order, degree, result = best

        print(
            f"  order={order} degree={degree}"
        )

        failures = verify_operator(
            B_ODD,
            D,
            result["solution"],
            order,
            degree,
        )

        boundary_checks = boundary_zero_residuals(
            B_ODD,
            D,
            result["solution"],
            order,
            degree,
        )

        boundary_ok = all_boundary_constraints_zero(
            boundary_checks
        )

        print_solution_blocks(
            result["solution"],
            order,
            degree,
        )

        print(
            f"  pointwise_exact={len(failures) == 0}"
        )
        print(
            f"  boundary_constraints_exact={boundary_ok}"
        )

        best_ok = (
            len(failures) == 0
            and boundary_ok
        )

    print()
    print("=" * 78)
    print("5. BOUNDARY FACTOR AUDIT OF BEST MODEL")
    print("=" * 78)

    if best is not None:
        order, degree, result = best

        boundary_factor_audit(
            B_ODD,
            D,
            result["solution"],
            order,
            degree,
        )
    else:
        print("  skipped: no unique exact model")

    cross_sector_test()

    print()
    print("=" * 78)
    print("6. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  Experiment 123 found that a generic degree<=2 p/r operator was
  not unique.

  Experiment 124 changes coordinates to

      d = D(p) - r,

  where D(p) is the exact terminal support of the residual row.

  This is the natural boundary-distance variable.

  The search therefore tests

      q_{p+2}(r)
        = sum_{a=0}^m F_a(p,d) q_p(r+a),

  with F_a polynomial in (p,d).

  The crucial distinction is:

      generic p/r law
          versus
      boundary-aware p/d law.

  A UNIQUE-EXACT low-order model with exact boundary zeros would be
  substantially stronger evidence for a genuine finite-support
  construction.

  A model that exists only after introducing many free parameters,
  or that fails the explicit boundary-zero constraints, is not
  promoted as structural.

  Everything is exact over Fraction arithmetic.
  No floating point.
  No SymPy.
  No extrapolation.
  """
    )

    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    print(f"  data_exact={data_ok}")
    print(f"  unique_boundary_model_found={best is not None}")
    print(f"  best_model_pointwise_exact={best_ok}")
    print(
        "  ALL BASIC CHECKS PASS="
        f"{data_ok and best is not None and best_ok}"
    )

    print()
    if data_ok and best is not None and best_ok:
        print("EXPERIMENT 124 COMPLETE")
    else:
        print("EXPERIMENT 124 COMPLETE — NO QUALIFIED UNIQUE BOUNDARY LAW")


if __name__ == "__main__":
    main()

