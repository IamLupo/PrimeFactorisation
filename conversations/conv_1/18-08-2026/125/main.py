from fractions import Fraction
from math import gcd


# ============================================================================
# EXPERIMENT 125 — EXACT BOUNDARY-DEFECT SHIFT-OPERATOR AUDIT
# ============================================================================
#
# Previous result:
#
#   q_{p+2}(r) = sum_a F_a(p,d) q_p(r+a)
#
# has no consistent low-degree boundary-aware solution, even with
#
#   d = D(p) - r.
#
# This experiment asks whether the obstruction is localized at the
# finite-support boundary.
#
# We therefore test
#
#   q_{p+2}(r)
#       =
#       sum_{a=0}^m F_a(p,d) q_p(r+a)
#       +
#       sum_{t=0}^{T} H_t(p) * I[d=t].
#
# Equivalently, the recurrence is homogeneous in the interior and
# receives explicit corrections only at the last T+1 boundary layers.
#
# The correction functions H_t(p) are deliberately restricted to
# low-degree polynomials in p.
#
# A positive result means:
#
#   interior propagation + finite boundary source
#
# is an exact construction mechanism.
#
# A negative result means even a boundary-localized defect does not
# explain the failure of the homogeneous operator at this complexity.
#
# Arithmetic:
#   exact Fraction only
#   no floating point
#   no SymPy
#   no extrapolation
#
# ============================================================================


# ----------------------------------------------------------------------------
# 1. EXACT B-ODD SECOND-LAYER DATA
# ----------------------------------------------------------------------------

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


# ----------------------------------------------------------------------------
# 2. EXACT LINEAR ALGEBRA
# ----------------------------------------------------------------------------

def matrix_rank(A):
    if not A:
        return 0

    M = [list(map(Fraction, row)) for row in A]
    rows = len(M)
    cols = len(M[0])

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

            f = M[r][col]
            if f == 0:
                continue

            for j in range(col, cols):
                M[r][j] -= f * M[rank][j]

        rank += 1
        col += 1

    return rank


def rref(A):
    if not A:
        return [], []

    M = [list(map(Fraction, row)) for row in A]
    rows = len(M)
    cols = len(M[0])

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

            f = M[r][col]
            if f == 0:
                continue

            for j in range(col, cols):
                M[r][j] -= f * M[rank][j]

        pivots.append(col)
        rank += 1
        col += 1

    return M, pivots


def solve_exact(A, b):
    if not A:
        return []

    aug = [
        list(map(Fraction, A[i])) + [Fraction(b[i])]
        for i in range(len(A))
    ]

    R, pivots = rref(aug)

    n = len(A[0])

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
# 3. BASIS FUNCTIONS
# ----------------------------------------------------------------------------

def basis_exponents(degree):
    exps = []

    for i in range(degree + 1):
        for j in range(degree + 1 - i):
            exps.append((i, j))

    return exps


def eval_pd_basis(p, d, degree):
    return [
        Fraction(p ** i) * Fraction(d ** j)
        for i, j in basis_exponents(degree)
    ]


def eval_p_basis(p, degree):
    return [
        Fraction(p ** i)
        for i in range(degree + 1)
    ]


def get_q(row, r):
    if 0 <= r < len(row):
        return row[r]
    return Fraction(0)


# ----------------------------------------------------------------------------
# 4. UNKNOWN LAYOUT
# ----------------------------------------------------------------------------
#
# Main operator:
#
#   F_a(p,d), a=0..order
#
# Boundary defects:
#
#   H_t(p), t=0..T
#
# Each H_t is activated only when d=t.
#
# ----------------------------------------------------------------------------

def operator_unknown_layout(order, operator_degree, boundary_T,
                            boundary_degree):
    layout = []

    # F_a blocks
    pd_basis = basis_exponents(operator_degree)

    for a in range(order + 1):
        for basis_idx in range(len(pd_basis)):
            layout.append(("F", a, basis_idx))

    # H_t blocks
    for t in range(boundary_T + 1):
        for j in range(boundary_degree + 1):
            layout.append(("H", t, j))

    return layout


def evaluate_unknown_vector(
    solution,
    layout,
    p,
    d,
    operator_degree,
    boundary_T,
    boundary_degree,
):
    total = Fraction(0)

    pd_basis = eval_pd_basis(p, d, operator_degree)
    p_basis = eval_p_basis(p, boundary_degree)

    for value, descriptor in zip(solution, layout):

        kind = descriptor[0]

        if kind == "F":
            _, a, idx = descriptor
            total += value * pd_basis[idx]

        else:
            _, t, idx = descriptor

            if d == t:
                total += value * p_basis[idx]

    return total


# ----------------------------------------------------------------------------
# 5. BUILD EXACT SYSTEM
# ----------------------------------------------------------------------------

def build_system(
    rows,
    Dmap,
    order,
    operator_degree,
    boundary_T,
    boundary_degree,
    step=2,
):
    layout = operator_unknown_layout(
        order,
        operator_degree,
        boundary_T,
        boundary_degree,
    )

    nvars = len(layout)

    equations = []
    rhs = []

    ps = sorted(rows)

    transitions = [
        (p, p + step)
        for p in ps
        if p + step in rows
    ]

    for p, p2 in transitions:

        Dp = Dmap[p]
        Dp2 = Dmap[p2]

        # Need enough r-values to test both ordinary propagation
        # and vanishing outside the source support.
        max_r = max(
            Dp + order + 2,
            Dp2 + 2,
        )

        for r in range(max_r + 1):

            d = Dp - r

            row_eq = [Fraction(0) for _ in range(nvars)]

            # F_a(p,d) * q_p(r+a)
            pd_basis = eval_pd_basis(
                p,
                d,
                operator_degree,
            )

            pd_basis_count = len(pd_basis)

            cursor = 0

            for a in range(order + 1):

                qvalue = get_q(rows[p], r + a)

                if qvalue != 0:
                    for j, bval in enumerate(pd_basis):
                        row_eq[cursor + j] += qvalue * bval

                cursor += pd_basis_count

            # H_t boundary terms.
            p_basis = eval_p_basis(
                p,
                boundary_degree,
            )

            for t in range(boundary_T + 1):
                if d != t:
                    continue

                for j, bval in enumerate(p_basis):
                    row_eq[cursor + t * len(p_basis) + j] += bval

            target = get_q(rows[p2], r)

            equations.append(row_eq)
            rhs.append(target)

    return equations, rhs, layout


# ----------------------------------------------------------------------------
# 6. OPERATOR SEARCH
# ----------------------------------------------------------------------------

def search_model(
    rows,
    Dmap,
    order,
    operator_degree,
    boundary_T,
    boundary_degree,
):
    A, b, layout = build_system(
        rows,
        Dmap,
        order,
        operator_degree,
        boundary_T,
        boundary_degree,
    )

    nvars = len(layout)

    rank_A = matrix_rank(A)

    augmented = [
        A[i] + [b[i]]
        for i in range(len(A))
    ]

    rank_aug = matrix_rank(augmented)

    if rank_A != rank_aug:
        return {
            "consistent": False,
            "unique": False,
            "rank": rank_A,
            "unknowns": nvars,
            "nullity": None,
            "solution": None,
            "layout": layout,
        }

    solution = solve_exact(A, b)

    if solution is None:
        return {
            "consistent": False,
            "unique": False,
            "rank": rank_A,
            "unknowns": nvars,
            "nullity": None,
            "solution": None,
            "layout": layout,
        }

    nullity = nvars - rank_A
    unique = nullity == 0

    return {
        "consistent": True,
        "unique": unique,
        "rank": rank_A,
        "unknowns": nvars,
        "nullity": nullity,
        "solution": solution,
        "layout": layout,
    }


# ----------------------------------------------------------------------------
# 7. POINTWISE VERIFICATION
# ----------------------------------------------------------------------------

def verify_model(
    rows,
    Dmap,
    result,
    order,
    operator_degree,
    boundary_T,
    boundary_degree,
):
    solution = result["solution"]
    layout = result["layout"]

    failures = []

    for p in sorted(rows):

        p2 = p + 2
        if p2 not in rows:
            continue

        Dp2 = Dmap[p2]

        # Verify a broad finite range, including points outside the
        # target support where the target is exactly zero.
        for r in range(Dp2 + order + boundary_T + 3):

            d = Dmap[p] - r

            rhs = Fraction(0)

            pd_basis = eval_pd_basis(
                p,
                d,
                operator_degree,
            )

            cursor = 0

            for a in range(order + 1):

                qvalue = get_q(rows[p], r + a)

                if qvalue != 0:
                    block = solution[
                        cursor:cursor + len(pd_basis)
                    ]

                    rhs += qvalue * sum(
                        block[j] * pd_basis[j]
                        for j in range(len(pd_basis))
                    )

                cursor += len(pd_basis)

            p_basis = eval_p_basis(
                p,
                boundary_degree,
            )

            for t in range(boundary_T + 1):

                if d != t:
                    continue

                start = cursor + t * len(p_basis)
                block = solution[
                    start:start + len(p_basis)
                ]

                rhs += sum(
                    block[j] * p_basis[j]
                    for j in range(len(p_basis))
                )

            lhs = get_q(rows[p2], r)

            if lhs != rhs:
                failures.append(
                    (p, p2, r, lhs, rhs)
                )

    return failures


# ----------------------------------------------------------------------------
# 8. EXTRACT / PRINT COMPONENTS
# ----------------------------------------------------------------------------

def extract_F_blocks(
    result,
    order,
    operator_degree,
):
    solution = result["solution"]
    layout = result["layout"]

    size = len(basis_exponents(operator_degree))

    F = []

    pos = 0
    for _a in range(order + 1):
        F.append(solution[pos:pos + size])
        pos += size

    return F


def extract_H_blocks(
    result,
    order,
    operator_degree,
    boundary_T,
    boundary_degree,
):
    solution = result["solution"]
    layout = result["layout"]

    Fsize = len(basis_exponents(operator_degree))
    Hsize = boundary_degree + 1

    start = (order + 1) * Fsize

    H = []

    for t in range(boundary_T + 1):
        H.append(
            solution[
                start + t * Hsize:
                start + (t + 1) * Hsize
            ]
        )

    return H


def primitive_signature(coeffs):
    if not coeffs:
        return []

    den_lcm = 1

    for c in coeffs:
        den_lcm = den_lcm * c.denominator // gcd(
            den_lcm,
            c.denominator,
        )

    ints = [
        int(c * den_lcm)
        for c in coeffs
    ]

    g = 0
    for n in ints:
        g = gcd(g, abs(n))

    if g == 0:
        return [0] * len(ints)

    ints = [n // g for n in ints]

    for n in ints:
        if n != 0:
            if n < 0:
                ints = [-x for x in ints]
            break

    return ints


def coefficient_digit_complexity(coeffs):
    if not coeffs:
        return (0, 0)

    max_num = max(
        len(str(abs(c.numerator)))
        for c in coeffs
    )

    max_den = max(
        len(str(c.denominator))
        for c in coeffs
    )

    return max_num, max_den


def print_blocks(
    result,
    order,
    operator_degree,
    boundary_T,
    boundary_degree,
):
    F = extract_F_blocks(
        result,
        order,
        operator_degree,
    )

    H = extract_H_blocks(
        result,
        order,
        operator_degree,
        boundary_T,
        boundary_degree,
    )

    pd_basis = basis_exponents(operator_degree)

    print("  OPERATOR BLOCKS")

    for a, block in enumerate(F):

        terms = []

        for coeff, (i, j) in zip(block, pd_basis):
            if coeff == 0:
                continue

            if i == 0 and j == 0:
                mon = "1"
            elif i == 0:
                mon = "d" if j == 1 else f"d^{j}"
            elif j == 0:
                mon = "p" if i == 1 else f"p^{i}"
            else:
                pm = "p" if i == 1 else f"p^{i}"
                dm = "d" if j == 1 else f"d^{j}"
                mon = f"{pm}*{dm}"

            terms.append(f"({coeff})*{mon}")

        expr = " + ".join(terms) if terms else "0"

        print(
            f"    F_{a}(p,d) = {expr}"
        )
        print(
            f"      primitive_signature="
            f"{primitive_signature(block)}"
        )
        print(
            f"      digit_complexity="
            f"{coefficient_digit_complexity(block)}"
        )

    print("  BOUNDARY DEFECT BLOCKS")

    for t, block in enumerate(H):

        terms = []

        for j, coeff in enumerate(block):

            if coeff == 0:
                continue

            mon = "1" if j == 0 else (
                "p" if j == 1 else f"p^{j}"
            )

            terms.append(
                f"({coeff})*{mon}"
            )

        expr = " + ".join(terms) if terms else "0"

        print(
            f"    H_{t}(p) = {expr}"
        )
        print(
            f"      primitive_signature="
            f"{primitive_signature(block)}"
        )
        print(
            f"      digit_complexity="
            f"{coefficient_digit_complexity(block)}"
        )


# ----------------------------------------------------------------------------
# 9. HOMOGENEOUS-VS-DEFECT COMPARISON
# ----------------------------------------------------------------------------

def homogeneous_equivalent_search(rows, Dmap):
    best = None

    for order in (1, 2, 3):
        for degree in (0, 1, 2, 3):

            result = search_model(
                rows,
                Dmap,
                order,
                degree,
                boundary_T=-1,
                boundary_degree=0,
            )

            # boundary_T=-1 is not used by the system builder directly,
            # so this function is only a descriptive placeholder.
            del result

    return best


# ----------------------------------------------------------------------------
# 10. MAIN SEARCH
# ----------------------------------------------------------------------------

def systematic_search(rows, Dmap):

    records = []

    print("=" * 78)
    print("EXPERIMENT 125 — EXACT BOUNDARY-DEFECT SHIFT-OPERATOR AUDIT")
    print("=" * 78)

    print()
    print("1. EXACT DATA")
    print("=" * 78)

    for p in sorted(rows):
        print(
            f"  p={p}: degree={len(rows[p]) - 1} "
            f"D={Dmap[p]} "
            f"entries={len(rows[p])}"
        )

    print()
    print("2. SYSTEMATIC BOUNDARY-DEFECT SEARCH")
    print("=" * 78)

    # The search is deliberately small.
    #
    # order        = 1,2,3
    # F degree     = 0,1,2
    # boundary T   = 0,1,2
    # H degree     = 0,1
    #
    for order in (1, 2, 3):

        for op_degree in (0, 1, 2):

            for boundary_T in (0, 1, 2):

                for boundary_degree in (0, 1):

                    result = search_model(
                        rows,
                        Dmap,
                        order,
                        op_degree,
                        boundary_T,
                        boundary_degree,
                    )

                    records.append(
                        (
                            order,
                            op_degree,
                            boundary_T,
                            boundary_degree,
                            result,
                        )
                    )

                    print(
                        f"  order={order} "
                        f"op_degree={op_degree} "
                        f"boundary_T={boundary_T} "
                        f"H_degree={boundary_degree}: "
                        f"consistent={result['consistent']} "
                        f"rank={result['rank']} "
                        f"unknowns={result['unknowns']} "
                        f"nullity={result['nullity']}"
                    )

    return records


# ----------------------------------------------------------------------------
# 11. SELECT MINIMAL QUALIFIED MODEL
# ----------------------------------------------------------------------------

def select_minimal(records, rows, Dmap):

    qualified = []

    for (
        order,
        op_degree,
        boundary_T,
        boundary_degree,
        result,
    ) in records:

        if not (
            result["consistent"]
            and result["unique"]
            and result["solution"] is not None
        ):
            continue

        failures = verify_model(
            rows,
            Dmap,
            result,
            order,
            op_degree,
            boundary_T,
            boundary_degree,
        )

        if failures:
            continue

        # Complexity ordering:
        #   first operator order,
        #   then operator degree,
        #   then number of boundary layers,
        #   then boundary degree.
        score = (
            order,
            op_degree,
            boundary_T,
            boundary_degree,
        )

        qualified.append(
            (
                score,
                order,
                op_degree,
                boundary_T,
                boundary_degree,
                result,
            )
        )

    if not qualified:
        return None

    return min(
        qualified,
        key=lambda item: item[0],
    )


# ----------------------------------------------------------------------------
# 12. BOUNDARY-DEFECT SUPPORT AUDIT
# ----------------------------------------------------------------------------

def boundary_support_audit(
    rows,
    Dmap,
    result,
    order,
    op_degree,
    boundary_T,
    boundary_degree,
):
    F = extract_F_blocks(
        result,
        order,
        op_degree,
    )

    H = extract_H_blocks(
        result,
        order,
        op_degree,
        boundary_T,
        boundary_degree,
    )

    print()
    print("=" * 78)
    print("4. EXACT BOUNDARY-DEFECT SUPPORT AUDIT")
    print("=" * 78)

    for t, block in enumerate(H):

        all_zero = all(c == 0 for c in block)

        print(
            f"  H_{t}: all_zero={all_zero} "
            f"coefficients={block}"
        )


# ----------------------------------------------------------------------------
# 13. INTERIOR-ONLY AUDIT
# ----------------------------------------------------------------------------

def interior_residual_audit(
    rows,
    Dmap,
    result,
    order,
    op_degree,
    boundary_T,
    boundary_degree,
):
    """
    Recheck the operator after deleting all points with d <= boundary_T.

    This determines whether the operator block F_a itself is exact
    in the genuine interior, with all mismatch confined to the boundary.
    """
    solution = result["solution"]

    F = extract_F_blocks(
        result,
        order,
        op_degree,
    )

    failures = []

    for p in sorted(rows):

        p2 = p + 2
        if p2 not in rows:
            continue

        Dp2 = Dmap[p2]

        for r in range(Dp2 + 1):

            d = Dmap[p] - r

            if d <= boundary_T:
                continue

            lhs = get_q(rows[p2], r)

            rhs = Fraction(0)

            pd_basis = eval_pd_basis(
                p,
                d,
                op_degree,
            )

            for a, block in enumerate(F):
                qvalue = get_q(
                    rows[p],
                    r + a,
                )

                if qvalue == 0:
                    continue

                value = sum(
                    block[j] * pd_basis[j]
                    for j in range(len(pd_basis))
                )

                rhs += qvalue * value

            if lhs != rhs:
                failures.append(
                    (p, p2, r, d, lhs, rhs)
                )

    return failures


# ----------------------------------------------------------------------------
# 14. MAIN
# ----------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 125 — EXACT BOUNDARY-DEFECT SHIFT-OPERATOR AUDIT")
    print("=" * 78)

    data_ok = True

    for p in sorted(B_ODD):

        degree = len(B_ODD[p]) - 1

        if degree != D[p]:
            data_ok = False

        if any(not isinstance(x, Fraction) for x in B_ODD[p]):
            data_ok = False

    print()
    print("=" * 78)
    print("1. DATA VALIDATION")
    print("=" * 78)

    for p in sorted(B_ODD):
        print(
            f"  p={p}: "
            f"D={D[p]} "
            f"degree={len(B_ODD[p]) - 1} "
            f"exact={len(B_ODD[p]) - 1 == D[p]}"
        )

    print(
        f"  data_exact={data_ok}"
    )

    records = systematic_search(
        B_ODD,
        D,
    )

    best = select_minimal(
        records,
        B_ODD,
        D,
    )

    print()
    print("=" * 78)
    print("3. LOWEST-COMPLEXITY UNIQUE EXACT MODEL")
    print("=" * 78)

    best_exact = False

    if best is None:

        print("  NONE")

    else:

        (
            score,
            order,
            op_degree,
            boundary_T,
            boundary_degree,
            result,
        ) = best

        print(
            f"  selected score={score}"
        )
        print(
            f"  order={order}"
        )
        print(
            f"  operator_degree={op_degree}"
        )
        print(
            f"  boundary_T={boundary_T}"
        )
        print(
            f"  boundary_degree={boundary_degree}"
        )

        failures = verify_model(
            B_ODD,
            D,
            result,
            order,
            op_degree,
            boundary_T,
            boundary_degree,
        )

        print(
            f"  pointwise_exact={len(failures) == 0}"
        )

        print_blocks(
            result,
            order,
            op_degree,
            boundary_T,
            boundary_degree,
        )

        interior_failures = interior_residual_audit(
            B_ODD,
            D,
            result,
            order,
            op_degree,
            boundary_T,
            boundary_degree,
        )

        print(
            f"  interior_exact="
            f"{len(interior_failures) == 0}"
        )

        best_exact = (
            len(failures) == 0
            and len(interior_failures) == 0
        )

    if best is not None:

        (
            _score,
            order,
            op_degree,
            boundary_T,
            boundary_degree,
            result,
        ) = best

        boundary_support_audit(
            B_ODD,
            D,
            result,
            order,
            op_degree,
            boundary_T,
            boundary_degree,
        )

    print()
    print("=" * 78)
    print("5. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  Experiment 124 showed that no homogeneous boundary-aware
  p/d operator exists at the tested low complexity.

  Experiment 125 therefore introduces the smallest natural
  inhomogeneous extension:

      q_{p+2}(r)
        =
        sum_a F_a(p,d) q_p(r+a)
        +
        sum_t H_t(p) [d=t],

  where

      d = D(p)-r.

  The term [d=t] is an exact boundary indicator.

  This separates two possible mechanisms:

      INTERIOR:
          reusable propagation operator F_a;

      BOUNDARY:
          finite source terms H_t.

  A UNIQUE-EXACT model with a low-order F and a small number
  of nonzero H_t would be a substantially stronger structural
  signal than an unrestricted interpolation.

  The interior audit is important: it verifies that the recurrence
  genuinely propagates through all non-boundary points rather than
  merely fitting the finite dataset through the defect variables.

  All arithmetic is exact Fraction arithmetic.
  No floating point.
  No SymPy.
  No extrapolation.
  """
    )

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  data_exact={data_ok}"
    )
    print(
        f"  unique_boundary_defect_model="
        f"{best is not None}"
    )
    print(
        f"  selected_model_exact={best_exact}"
    )
    print(
        "  ALL BASIC CHECKS PASS="
        f"{data_ok and best is not None and best_exact}"
    )

    print()

    if data_ok and best is not None and best_exact:
        print("EXPERIMENT 125 COMPLETE")
    else:
        print(
            "EXPERIMENT 125 COMPLETE — "
            "NO QUALIFIED UNIQUE BOUNDARY-DEFECT LAW"
        )


if __name__ == "__main__":
    main()

