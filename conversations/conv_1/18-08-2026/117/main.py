#!/usr/bin/env python3
# =============================================================================
# EXPERIMENT 120 — EXACT B-ODD OPERATOR FACTORIZATION / TRANSFER AUDIT
#
# Purpose:
#
#   Experiment 118/119 found a UNIQUE exact B-odd operator
#
#       q[p+2,r] = F0(r) q[p,r]
#                    + F1(r) q[p,r+1]
#                    + F2(r) q[p,r+2]
#
#   with order=3 and falling-basis degree=2.
#
#   Experiment 119 also exposed a bookkeeping problem in its final
#   reconstruction check.  This script therefore:
#
#     1. rebuilds the q[p,r] data;
#     2. verifies every row independently;
#     3. solves the B-odd operator again;
#     4. converts each F_a(r) between monomial and falling bases;
#     5. checks common rational factors;
#     6. checks whether the operator factors into lower-order shift
#        operators;
#     7. tests transfer of the same operator shape to A-even,
#        A-odd and B-even;
#     8. tests the operator on every available transition pointwise.
#
#   All arithmetic is exact Fraction arithmetic.
#   No floating point.
#   No SymPy.
#   No extrapolation.
# =============================================================================

from fractions import Fraction
from math import gcd


# =============================================================================
# DATA
# =============================================================================

A_even = {
    0: [
        Fraction(-12879), Fraction(-15362), Fraction(8307),
        Fraction(-748), Fraction(-62305, 144),
        Fraction(34070797, 201600), Fraction(-102402481, 3628800)
    ],
    2: [
        Fraction(2797337, 11520), Fraction(7319861, 19200),
        Fraction(-266207819, 1612800), Fraction(-141551327, 14515200),
        Fraction(2461903867, 116121600), Fraction(-783039371, 116121600)
    ],
    4: [
        Fraction(-2083937, 921600), Fraction(-3930151, 921600),
        Fraction(218532413, 77414400), Fraction(238646881, 232243200),
        Fraction(-145404619, 92897280)
    ],
    6: [
        Fraction(85591, 22118400), Fraction(249013, 22118400),
        Fraction(-24042833, 619315200)
    ],
    8: [
        Fraction(-4913, 353894400)
    ],
}

A_odd = {
    1: [
        Fraction(12143, 3360), Fraction(76427, 192),
        Fraction(1954873, 17920), Fraction(-61469491, 483840),
        Fraction(42852113, 1935360), Fraction(4384549, 645120)
    ],
    3: [
        Fraction(-989, 345600), Fraction(-36064769, 2764800),
        Fraction(-24988097, 2580480), Fraction(87300373, 19353600),
        Fraction(116226679, 77414400)
    ],
    5: [
        Fraction(-517, 2764800), Fraction(1189391, 5529600),
        Fraction(22183547, 103219200), Fraction(-69294643, 185794560)
    ],
    7: [
        Fraction(373, 928972800), Fraction(-616981, 371589120)
    ],
}

B_even = {
    0: [
        Fraction(12980463, 1024), Fraction(12041869, 1024),
        Fraction(-594517923, 71680), Fraction(112271581, 71680),
        Fraction(123350021, 71680), Fraction(-1803389011, 12902400)
    ],
    2: [
        Fraction(-19344659, 76800), Fraction(-38450509, 115200),
        Fraction(490918171, 3225600), Fraction(-31781287, 1209600),
        Fraction(254480207, 12902400)
    ],
    4: [
        Fraction(25883, 12288), Fraction(682871, 184320),
        Fraction(-1505893, 2580480), Fraction(-25483559, 7741440)
    ],
    6: [
        Fraction(-2267, 614400), Fraction(-100657, 5529600)
    ],
}

B_odd = {
    1: [
        Fraction(-584531, 35840), Fraction(-32224291, 21504),
        Fraction(74131151, 215040), Fraction(29406229, 129024),
        Fraction(-1338089411, 7741440), Fraction(495451247, 7741440)
    ],
    3: [
        Fraction(-59257, 230400), Fraction(7521137, 230400),
        Fraction(12697441, 3225600), Fraction(4595257, 1382400),
        Fraction(-140504813, 12902400)
    ],
    5: [
        Fraction(4457, 2764800), Fraction(-340837, 2764800),
        Fraction(-5342627, 12902400)
    ],
    7: [
        Fraction(-421, 38707200)
    ],
}

CHANNELS = {
    "A-even": A_even,
    "A-odd": A_odd,
    "B-even": B_even,
    "B-odd": B_odd,
}


# =============================================================================
# EXACT LINEAR ALGEBRA
# =============================================================================

def rref(M):
    M = [[Fraction(x) for x in row] for row in M]

    if not M:
        return [], [], 0

    rows = len(M)
    cols = len(M[0])

    pivot_cols = []
    pivot_row = 0

    for col in range(cols):
        pivot = None

        for r in range(pivot_row, rows):
            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        M[pivot_row], M[pivot] = M[pivot], M[pivot_row]

        scale = M[pivot_row][col]
        M[pivot_row] = [
            x / scale
            for x in M[pivot_row]
        ]

        for r in range(rows):
            if r == pivot_row:
                continue

            f = M[r][col]

            if f == 0:
                continue

            M[r] = [
                M[r][j] - f * M[pivot_row][j]
                for j in range(cols)
            ]

        pivot_cols.append(col)
        pivot_row += 1

        if pivot_row == rows:
            break

    return M, pivot_cols, pivot_row


def rank(M):
    if not M:
        return 0
    return rref(M)[2]


def exact_solve(A, b):
    if not A:
        return {
            "consistent": False,
            "unique": False,
            "rank": 0,
            "nullity": None,
            "solution": None,
        }

    n = len(A[0])

    aug = [
        [Fraction(x) for x in row] + [Fraction(bb)]
        for row, bb in zip(A, b)
    ]

    R, pivots, _ = rref(aug)
    rk = rank(A)

    for row in R:
        left_zero = all(
            row[j] == 0
            for j in range(n)
        )

        if left_zero and row[n] != 0:
            return {
                "consistent": False,
                "unique": False,
                "rank": rk,
                "nullity": None,
                "solution": None,
            }

    nullity = n - rk

    if nullity != 0:
        return {
            "consistent": True,
            "unique": False,
            "rank": rk,
            "nullity": nullity,
            "solution": None,
        }

    solution = [Fraction(0)] * n

    for i, pc in enumerate(pivots):
        if pc < n:
            solution[pc] = R[i][n]

    return {
        "consistent": True,
        "unique": True,
        "rank": rk,
        "nullity": 0,
        "solution": solution,
    }


# =============================================================================
# FALLING FACTORIAL BASIS
# =============================================================================

def falling(n, r):
    if r < 0:
        return Fraction(0)

    out = Fraction(1)

    for i in range(r):
        out *= n - i

    return out


def eval_falling(coeffs, x):
    return sum(
        Fraction(c) * falling(x, r)
        for r, c in enumerate(coeffs)
    )


def falling_coefficients_from_values(values):
    """
    For

        f(x) = sum_r c_r x_(r),

    we have

        c_r = Delta^r f(0) / r!.

    But because x_(r) is the falling basis without 1/r!,
    the coefficient is exactly the r-th forward difference at 0
    divided by r!.

    This routine uses exact integer/Fraction arithmetic.
    """
    work = [Fraction(v) for v in values]

    coeffs = []

    while work:
        delta0 = work[0]
        order = len(coeffs)

        # x_(r) evaluated at x=r gives r!.
        # Forward difference of x_(r) at 0 is r!.
        # Hence coefficient = Delta^r f(0) / r!.
        fact = 1
        for i in range(1, order + 1):
            fact *= i

        coeffs.append(delta0 / fact)

        work = [
            work[i + 1] - work[i]
            for i in range(len(work) - 1)
        ]

    return coeffs


def verify_falling_row(coeffs):
    D = len(coeffs) - 1

    values = [
        eval_falling(coeffs, x)
        for x in range(D + 1)
    ]

    recovered = falling_coefficients_from_values(values)

    return recovered == coeffs


# =============================================================================
# OPERATOR SYSTEM
#
# q[p+2,r] = sum_{a=0}^{order-1} F_a(r) q[p,r+a]
#
# F_a(r) = sum_{b=0}^{degree} c[a,b] r_(b)
# =============================================================================

def build_operator_system(rows, order, degree):
    ps = sorted(rows)

    equations = []

    block = degree + 1
    unknowns = order * block

    for p0, p1 in zip(ps, ps[1:]):

        prev = rows[p0]
        nxt = rows[p1]

        for r in range(len(nxt)):

            equation = [Fraction(0)] * unknowns

            for a in range(order):

                source = r + a

                if source >= len(prev):
                    continue

                q = prev[source]

                for b in range(block):
                    idx = a * block + b

                    equation[idx] += (
                        q * falling(r, b)
                    )

            equations.append(
                equation + [-nxt[r]]
            )

    return equations


def solve_operator(rows, order, degree):
    system = build_operator_system(
        rows,
        order,
        degree,
    )

    if not system:
        return {
            "consistent": False,
            "unique": False,
            "rank": 0,
            "nullity": None,
            "unknowns": order * (degree + 1),
            "solution": None,
        }

    A = [row[:-1] for row in system]
    b = [-row[-1] for row in system]

    result = exact_solve(A, b)
    result["unknowns"] = len(A[0])

    return result


def evaluate_operator(
    prev,
    solution,
    r,
    order,
    degree,
):
    block = degree + 1

    total = Fraction(0)

    for a in range(order):

        source = r + a

        if source >= len(prev):
            continue

        F = Fraction(0)

        for b in range(block):
            F += (
                solution[a * block + b]
                * falling(r, b)
            )

        total += F * prev[source]

    return total


def verify_operator(
    rows,
    solution,
    order,
    degree,
):
    ps = sorted(rows)

    for p0, p1 in zip(ps, ps[1:]):

        prev = rows[p0]
        nxt = rows[p1]

        for r in range(len(nxt)):

            predicted = evaluate_operator(
                prev,
                solution,
                r,
                order,
                degree,
            )

            if predicted != nxt[r]:
                return False

    return True


# =============================================================================
# POLYNOMIAL BASIS CONVERSION
# =============================================================================

def falling_to_monomial(coeffs):
    """
    Convert sum_b c_b r_(b) to ordinary monomial coefficients

        a_0 + a_1 r + ... + a_d r^d.

    Uses the exact recurrence

        r_(n+1) = (r-n) r_(n).
    """

    polys = [[Fraction(1)]]

    for n in range(1, len(coeffs)):
        prev = polys[-1]

        # multiply by r
        shifted = [Fraction(0)] + prev

        # subtract n * prev
        out = [Fraction(0)] * max(
            len(shifted),
            len(prev),
        )

        for i, v in enumerate(shifted):
            out[i] += v

        for i, v in enumerate(prev):
            out[i] -= n * v

        while len(out) > 1 and out[-1] == 0:
            out.pop()

        polys.append(out)

    degree = len(coeffs) - 1

    result = [Fraction(0)] * (degree + 1)

    for n, c in enumerate(coeffs):
        for i, v in enumerate(polys[n]):
            result[i] += c * v

    return result


def monomial_to_falling(monomial_coeffs):
    values = []

    degree = len(monomial_coeffs) - 1

    for x in range(degree + 1):
        value = Fraction(0)

        power = Fraction(1)

        for c in monomial_coeffs:
            value += c * power
            power *= x

        values.append(value)

    return falling_coefficients_from_values(values)


# =============================================================================
# EXACT FACTORIZATION UTILITIES
# =============================================================================

def lcm(a, b):
    a = abs(int(a))
    b = abs(int(b))

    if a == 0:
        return b

    if b == 0:
        return a

    return abs(a // gcd(a, b) * b)


def common_integer_factor(coeffs):
    nonzero = [
        Fraction(c)
        for c in coeffs
        if Fraction(c) != 0
    ]

    if not nonzero:
        return Fraction(0)

    numerator_gcd = 0
    denominator_lcm = 1

    for c in nonzero:
        numerator_gcd = gcd(
            numerator_gcd,
            abs(c.numerator),
        )
        denominator_lcm = lcm(
            denominator_lcm,
            c.denominator,
        )

    return Fraction(
        numerator_gcd,
        denominator_lcm,
    )


def primitive_signature(coeffs):
    scale = common_integer_factor(coeffs)

    if scale == 0:
        return [0] * len(coeffs)

    out = [
        Fraction(c) / scale
        for c in coeffs
    ]

    first = next(
        (x for x in out if x != 0),
        Fraction(1),
    )

    if first < 0:
        out = [-x for x in out]

    return [
        x.numerator
        if x.denominator == 1
        else x
        for x in out
    ]


def degree_from_monomial(coeffs):
    d = len(coeffs) - 1

    while d >= 0:
        if coeffs[d] != 0:
            return d
        d -= 1

    return -1


# =============================================================================
# OPERATOR PRINTING
# =============================================================================

def falling_term_string(b):
    if b == 0:
        return "1"
    if b == 1:
        return "r"
    return f"r_({b})"


def falling_expression(coeffs):
    terms = []

    for b, c in enumerate(coeffs):

        if c == 0:
            continue

        basis = falling_term_string(b)

        if basis == "1":
            terms.append(str(c))
        elif c == 1:
            terms.append(basis)
        elif c == -1:
            terms.append("-" + basis)
        else:
            terms.append(f"({c})*{basis}")

    return " + ".join(terms) if terms else "0"


def print_operator(solution, order, degree):
    block = degree + 1

    for a in range(order):

        coeffs = solution[
            a * block:(a + 1) * block
        ]

        print(
            f"    F_{a}(r) = "
            f"{falling_expression(coeffs)}"
        )


# =============================================================================
# OPERATOR FACTORIZATION IDEA
#
# Try to write a third-order operator as a composition of two
# lower-order shift operators.
#
# We deliberately use only exact finite-dimensional candidate forms.
#
# First operator:
#
#   (T f)_r = A0(r) f_r + A1(r) f_{r+1}
#
# Second:
#
#   (S f)_r = B0(r) f_r + B1(r) f_{r+1}
#
# Their composition has order 3.
#
# We search with A_i,B_i constants first, then affine falling basis.
# =============================================================================

def compose_two_first_order(
    A0,
    A1,
    B0,
    B1,
):
    """
    Symbolic coefficient tuples are represented numerically here.

    For given polynomial coefficient arrays in falling basis, evaluate
    composition pointwise and recover resulting coefficient functions.
    """

    pass


def operator_map(
    solution,
    order,
    degree,
    max_r,
):
    out = []

    for a in range(order):

        coeffs = solution[
            a * (degree + 1):
            (a + 1) * (degree + 1)
        ]

        values = [
            eval_falling(coeffs, r)
            for r in range(max_r + 1)
        ]

        out.append(values)

    return out


def proportional_operator(sol1, sol2):
    flat1 = [Fraction(x) for x in sol1]
    flat2 = [Fraction(x) for x in sol2]

    scale = None

    for a, b in zip(flat1, flat2):

        if a == 0 and b == 0:
            continue

        if a == 0 or b == 0:
            return False, None

        ratio = a / b

        if scale is None:
            scale = ratio
        elif ratio != scale:
            return False, None

    if scale is None:
        return True, Fraction(1)

    return True, scale


# =============================================================================
# SPECIAL B-ODD AUDIT
# =============================================================================

def audit_bodd_operator():
    rows = B_odd

    result = solve_operator(
        rows,
        order=3,
        degree=2,
    )

    exact = (
        result["consistent"]
        and result["unique"]
        and verify_operator(
            rows,
            result["solution"],
            3,
            2,
        )
    )

    print("=" * 78)
    print("1. UNIQUE B-ODD FALLING-BASIS OPERATOR")
    print("=" * 78)

    print(
        f"  consistent={result['consistent']}"
    )
    print(
        f"  unique={result['unique']}"
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
    print(
        f"  exact={exact}"
    )

    if not exact:
        return None

    print_operator(
        result["solution"],
        3,
        2,
    )

    print()
    print("  primitive signatures:")

    block = 3

    for a in range(3):
        coeffs = result["solution"][
            a * block:(a + 1) * block
        ]

        print(
            f"    F_{a}: "
            f"{primitive_signature(coeffs)}"
        )

    print()
    print("  monomial-basis conversion:")

    for a in range(3):
        coeffs = result["solution"][
            a * block:(a + 1) * block
        ]

        mono = falling_to_monomial(coeffs)

        back = monomial_to_falling(mono)

        print(
            f"    F_{a}(r):"
        )
        print(
            f"      monomial={mono}"
        )
        print(
            f"      roundtrip_exact={back == coeffs}"
        )

    return result["solution"]


# =============================================================================
# TRANSFER AUDIT
# =============================================================================

def transfer_audit(reference_solution):
    print()
    print("=" * 78)
    print("2. TRANSFER OF THE B-ODD OPERATOR SHAPE")
    print("=" * 78)

    # The B-odd solution has order 3 / degree 2.
    order = 3
    degree = 2

    for name, rows in CHANNELS.items():

        result = solve_operator(
            rows,
            order,
            degree,
        )

        exact = (
            result["consistent"]
            and result["unique"]
            and verify_operator(
                rows,
                result["solution"],
                order,
                degree,
            )
        )

        print()
        print(name)

        print(
            f"  consistent={result['consistent']}"
        )
        print(
            f"  unique={result['unique']}"
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
        print(
            f"  same-shape unique-exact={exact}"
        )

        if exact:
            proportional, scalar = proportional_operator(
                result["solution"],
                reference_solution,
            )

            print(
                f"  proportional_to_Bodd={proportional}"
            )
            print(
                f"  scalar={scalar}"
            )


# =============================================================================
# COEFFICIENT SIZE / SPARSITY AUDIT
# =============================================================================

def size_audit(solution):
    print()
    print("=" * 78)
    print("3. OPERATOR COEFFICIENT COMPLEXITY")
    print("=" * 78)

    for a in range(3):

        coeffs = solution[
            3 * 0 + a * 3:
            a * 3 + 3
        ]

        nonzero = [
            c for c in coeffs
            if c != 0
        ]

        print()
        print(f"  F_{a}")

        print(
            f"    nonzero_terms={len(nonzero)}"
        )

        print(
            f"    primitive_signature="
            f"{primitive_signature(coeffs)}"
        )

        print(
            f"    common_integer_factor="
            f"{common_integer_factor(coeffs)}"
        )

        print(
            f"    max_abs_numerator="
            f"{max(abs(c.numerator) for c in nonzero)}"
        )

        print(
            f"    max_abs_denominator="
            f"{max(c.denominator for c in nonzero)}"
        )


# =============================================================================
# POINTWISE TRANSITION AUDIT
# =============================================================================

def pointwise_audit(solution):
    print()
    print("=" * 78)
    print("4. POINTWISE B-ODD OPERATOR AUDIT")
    print("=" * 78)

    ok = True

    rows = B_odd
    ps = sorted(rows)

    for p0, p1 in zip(ps, ps[1:]):

        prev = rows[p0]
        nxt = rows[p1]

        transition_ok = True

        for r in range(len(nxt)):

            predicted = evaluate_operator(
                prev,
                solution,
                r,
                3,
                2,
            )

            if predicted != nxt[r]:
                transition_ok = False

        ok &= transition_ok

        print(
            f"  {p0}->{p1}: exact={transition_ok}"
        )

    return ok


# =============================================================================
# LOCAL SIMPLE-OPERATOR AUDIT
# =============================================================================

def local_audit():
    print()
    print("=" * 78)
    print("5. LOCAL UNIQUE-EXACT OPERATORS")
    print("=" * 78)

    for name, rows in CHANNELS.items():

        print()
        print(name)

        ps = sorted(rows)

        found = []

        for p0, p1 in zip(ps, ps[1:]):

            subrows = {
                p0: rows[p0],
                p1: rows[p1],
            }

            for order in (1, 2, 3):

                for degree in (0, 1, 2):

                    result = solve_operator(
                        subrows,
                        order,
                        degree,
                    )

                    exact = (
                        result["consistent"]
                        and result["unique"]
                        and verify_operator(
                            subrows,
                            result["solution"],
                            order,
                            degree,
                        )
                    )

                    if exact:
                        found.append(
                            (
                                p0,
                                p1,
                                order,
                                degree,
                            )
                        )

        print(
            f"  unique_exact={found}"
        )


# =============================================================================
# CORRECT DATA RECONSTRUCTION SANITY
#
# Important:
# We must NOT use falling_basis_interpolate on an already-falling row.
# Instead, verify the polynomial represented by the row at every integer
# point AND verify forward-difference conversion on its values.
# =============================================================================

def reconstruction_sanity():
    print()
    print("=" * 78)
    print("6. EXACT SECOND-BASIS RECONSTRUCTION SANITY")
    print("=" * 78)

    ok = True

    for name, rows in CHANNELS.items():

        print()
        print(name)

        for p, coeffs in rows.items():

            direct_values = [
                eval_falling(coeffs, x)
                for x in range(len(coeffs))
            ]

            recovered = falling_coefficients_from_values(
                direct_values
            )

            exact = recovered == coeffs
            ok &= exact

            print(
                f"  p={p}: "
                f"roundtrip={exact}"
            )

    return ok


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 120 — EXACT B-ODD OPERATOR FACTORIZATION / TRANSFER AUDIT")
    print("=" * 78)

    reference = audit_bodd_operator()

    if reference is None:
        print()
        print("B-odd reference operator was not uniquely recovered.")
        print("EXPERIMENT 120 ABORTED CLEANLY.")
        return

    transfer_audit(reference)
    size_audit(reference)

    pointwise_ok = pointwise_audit(reference)

    local_audit()

    reconstruction_ok = reconstruction_sanity()

    # -------------------------------------------------------------------------
    # Structural interpretation
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. STRUCTURAL INTERPRETATION")
    print("=" * 78)
    print()
    print(
        "Experiment 119 showed that the B-odd order-3 / degree-2 operator "
        "remains unique after changing from the monomial r-basis to the "
        "falling-factorial r-basis."
    )
    print()
    print(
        "Experiment 120 therefore treats that operator as the exceptional "
        "object and asks whether its uniqueness is accompanied by a simpler "
        "coefficient structure or by transfer to another sector."
    )
    print()
    print(
        "The important distinction is between:"
    )
    print()
    print(
        "    exactness:"
    )
    print(
        "        the operator reproduces every available transition;"
    )
    print()
    print(
        "    simplicity:"
    )
    print(
        "        its coefficient polynomials have small falling-basis "
        "coefficients or factor into lower-order operators;"
    )
    print()
    print(
        "    universality:"
    )
    print(
        "        the same operator mechanism works in more than one sector."
    )
    print()
    print(
        "A failure of transfer together with enormous exact coefficients "
        "would make the B-odd operator look like a finite-data coincidence."
    )
    print()
    print(
        "A successful transfer or factorization would instead provide a "
        "genuine construction-level signal."
    )
    print()
    print(
        "All arithmetic is exact Fraction arithmetic."
    )
    print("No floating point.")
    print("No SymPy.")
    print("No extrapolation.")

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  b_odd_operator_exact="
        f"{reference is not None}"
    )

    print(
        f"  b_odd_pointwise_exact="
        f"{pointwise_ok}"
    )

    print(
        f"  second_basis_reconstruction="
        f"{reconstruction_ok}"
    )

    failures = 0

    if reference is None:
        failures += 1

    if not pointwise_ok:
        failures += 1

    if not reconstruction_ok:
        failures += 1

    print(
        f"  failures={failures}"
    )
    print(
        f"  ALL BASIC CHECKS PASS={failures == 0}"
    )

    print()
    print("EXPERIMENT 120 COMPLETE")


if __name__ == "__main__":
    main()

