#!/usr/bin/env python3

from fractions import Fraction
from math import gcd
from itertools import product


# =============================================================================
# EXPERIMENT 119
# EXACT FALLING-BASIS r-OPERATOR AUDIT
#
# Goal:
#
#   Experiment 118 found one unique exact operator:
#
#       q[p+2,r]
#         = F0(r) q[p,r]
#         + F1(r) q[p,r+1]
#         + F2(r) q[p,r+2]
#
#   for B-odd, with F_a represented as ordinary polynomials in r.
#
#   This experiment asks whether the SAME operator becomes much simpler
#   when each F_a(r) is represented in the falling-factorial basis:
#
#       F_a(r) = sum_b c[a,b] r_(b).
#
#   It also searches for exact operators in that basis across all sectors.
#
#   Arithmetic:
#       exact Fraction only
#       no floating point
#       no SymPy
#       no extrapolation
# =============================================================================


# -----------------------------------------------------------------------------
# DATA
# -----------------------------------------------------------------------------

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


# -----------------------------------------------------------------------------
# EXACT LINEAR ALGEBRA
# -----------------------------------------------------------------------------

def rref(matrix):
    if not matrix:
        return [], [], 0

    M = [[Fraction(x) for x in row] for row in matrix]
    nr = len(M)
    nc = len(M[0])

    pivots = []
    row = 0

    for col in range(nc):
        pivot = None

        for rr in range(row, nr):
            if M[rr][col] != 0:
                pivot = rr
                break

        if pivot is None:
            continue

        M[row], M[pivot] = M[pivot], M[row]

        pv = M[row][col]
        M[row] = [x / pv for x in M[row]]

        for rr in range(nr):
            if rr == row:
                continue

            factor = M[rr][col]

            if factor == 0:
                continue

            M[rr] = [
                M[rr][j] - factor * M[row][j]
                for j in range(nc)
            ]

        pivots.append(col)
        row += 1

        if row == nr:
            break

    return M, pivots, row


def matrix_rank(A):
    if not A:
        return 0
    return rref(A)[2]


def solve_exact(A, b):
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
    rk = matrix_rank(A)

    for row in R:
        if (
            all(row[j] == 0 for j in range(n))
            and row[n] != 0
        ):
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

    sol = [Fraction(0)] * n

    for i, pivot in enumerate(pivots):
        if pivot < n:
            sol[pivot] = R[i][n]

    return {
        "consistent": True,
        "unique": True,
        "rank": rk,
        "nullity": 0,
        "solution": sol,
    }


# -----------------------------------------------------------------------------
# FALLING FACTORIALS
# -----------------------------------------------------------------------------

def falling(n, r):
    if r < 0:
        return Fraction(0)

    out = Fraction(1)

    for i in range(r):
        out *= n - i

    return out


# -----------------------------------------------------------------------------
# BUILD FALLING-BASIS OPERATOR
#
#   q[p+2,r]
#     = sum_a F_a(r) q[p,r+a]
#
#   F_a(r) = sum_b c[a,b] r_(b)
# -----------------------------------------------------------------------------

def build_operator_system(
    rows,
    order,
    basis_degree,
):
    ps = sorted(rows)

    block = basis_degree + 1
    unknowns = order * block

    equations = []

    for p0, p1 in zip(ps, ps[1:]):
        prev = rows[p0]
        nxt = rows[p1]

        for r in range(len(nxt)):

            row = [Fraction(0)] * unknowns

            for a in range(order):

                source = r + a

                if source >= len(prev):
                    continue

                q = prev[source]

                for b in range(basis_degree + 1):
                    idx = a * block + b

                    row[idx] += (
                        q * falling(r, b)
                    )

            # operator - target = 0
            equations.append(row + [-nxt[r]])

    return equations


def solve_operator(rows, order, basis_degree):
    system = build_operator_system(
        rows,
        order,
        basis_degree,
    )

    if not system:
        return {
            "consistent": False,
            "unique": False,
            "rank": 0,
            "unknowns": order * (basis_degree + 1),
            "nullity": None,
            "solution": None,
        }

    A = [row[:-1] for row in system]
    b = [-row[-1] for row in system]

    result = solve_exact(A, b)
    result["unknowns"] = len(A[0])

    return result


def evaluate_operator(
    prev,
    solution,
    r,
    order,
    basis_degree,
):
    block = basis_degree + 1
    total = Fraction(0)

    for a in range(order):

        source = r + a

        if source >= len(prev):
            continue

        coeff = Fraction(0)
        base = a * block

        for b in range(basis_degree + 1):
            coeff += (
                solution[base + b]
                * falling(r, b)
            )

        total += coeff * prev[source]

    return total


def verify_operator(
    rows,
    order,
    basis_degree,
    solution,
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
                basis_degree,
            )

            if predicted != nxt[r]:
                return False

    return True


# -----------------------------------------------------------------------------
# RE-EXPRESS MONOMIAL POLYNOMIAL IN FALLING BASIS
# -----------------------------------------------------------------------------

def polynomial_values_from_falling(coeffs, max_r):
    out = []

    for r in range(max_r + 1):
        value = Fraction(0)

        for b, c in enumerate(coeffs):
            value += c * falling(r, b)

        out.append(value)

    return out


def falling_basis_interpolate(values):
    """
    Given values f(0),...,f(d), recover

        f(r) = sum_b c_b r_(b)

    exactly by forward differences.
    """
    d = len(values) - 1

    work = [Fraction(v) for v in values]
    coeffs = []

    while work:
        coeffs.append(work[0])

        work = [
            work[i + 1] - work[i]
            for i in range(len(work) - 1)
        ]

    return coeffs


# -----------------------------------------------------------------------------
# SIGNATURES / COMPLEXITY
# -----------------------------------------------------------------------------

def lcm(a, b):
    if a == 0 or b == 0:
        return 0
    return abs(a // gcd(a, b) * b)


def common_integer_factor(coeffs):
    nonzero = [Fraction(c) for c in coeffs if c != 0]

    if not nonzero:
        return Fraction(0)

    ng = 0
    dl = 1

    for c in nonzero:
        ng = gcd(ng, abs(c.numerator))
        dl = lcm(dl, c.denominator)

    return Fraction(ng, dl)


def primitive_signature(coeffs):
    nonzero = [Fraction(c) for c in coeffs if c != 0]

    if not nonzero:
        return [0] * len(coeffs)

    scale = common_integer_factor(nonzero)

    if scale == 0:
        return [0] * len(coeffs)

    out = [
        Fraction(c) / scale
        for c in coeffs
    ]

    # Make a deterministic sign choice.
    first = next(
        (x for x in out if x != 0),
        Fraction(1),
    )

    if first < 0:
        out = [-x for x in out]

    values = []

    for x in out:
        if x.denominator != 1:
            values.append(x)
        else:
            values.append(x.numerator)

    return values


def max_abs_numerator(coeffs):
    return max(
        abs(Fraction(c).numerator)
        for c in coeffs
        if c != 0
    ) if any(c != 0 for c in coeffs) else 0


def max_abs_denominator(coeffs):
    return max(
        Fraction(c).denominator
        for c in coeffs
        if c != 0
    ) if any(c != 0 for c in coeffs) else 1


# -----------------------------------------------------------------------------
# PRINTING
# -----------------------------------------------------------------------------

def operator_blocks(solution, order, degree):
    block = degree + 1
    return [
        solution[a * block:(a + 1) * block]
        for a in range(order)
    ]


def falling_expression(coeffs):
    terms = []

    for b, c in enumerate(coeffs):

        if c == 0:
            continue

        if b == 0:
            basis = "1"
        elif b == 1:
            basis = "r"
        else:
            basis = f"r_({b})"

        if basis == "1":
            term = str(c)
        elif c == 1:
            term = basis
        elif c == -1:
            term = "-" + basis
        else:
            term = f"({c})*{basis}"

        terms.append(term)

    return " + ".join(terms) if terms else "0"


def print_operator(solution, order, degree):
    blocks = operator_blocks(
        solution,
        order,
        degree,
    )

    for a, coeffs in enumerate(blocks):
        print(
            f"    F_{a}(r) = "
            f"{falling_expression(coeffs)}"
        )


# -----------------------------------------------------------------------------
# B-ODD MONOMIAL VS FALLING COMPARISON
#
# Rather than reproducing the enormous Experiment-118 coefficients,
# solve the same operator directly in both bases and compare exact
# coefficient complexity.
# -----------------------------------------------------------------------------

def compare_bodd():

    print("=" * 78)
    print("1. B-ODD MONOMIAL VS FALLING-BASIS OPERATOR")
    print("=" * 78)

    rows = B_odd

    candidates = []

    for degree in (0, 1, 2, 3):
        result = solve_operator(
            rows,
            order=3,
            basis_degree=degree,
        )

        exact = (
            result["unique"]
            and verify_operator(
                rows,
                3,
                degree,
                result["solution"],
            )
        )

        print(
            f"  falling-basis degree={degree}: "
            f"consistent={result['consistent']} "
            f"rank={result['rank']} "
            f"unknowns={result['unknowns']} "
            f"nullity={result['nullity']} "
            f"unique_exact={exact}"
        )

        if exact:
            solution = result["solution"]
            candidates.append(
                (degree, solution)
            )

            blocks = operator_blocks(
                solution,
                3,
                degree,
            )

            print(
                f"    max|numerator|="
                f"{max(max_abs_numerator(b) for b in blocks)}"
            )
            print(
                f"    max|denominator|="
                f"{max(max_abs_denominator(b) for b in blocks)}"
            )

            print_operator(
                solution,
                3,
                degree,
            )

    return candidates


# -----------------------------------------------------------------------------
# UNIVERSAL SEARCH
# -----------------------------------------------------------------------------

def universal_search():

    print()
    print("=" * 78)
    print("2. UNIVERSAL FALLING-BASIS r-OPERATOR SEARCH")
    print("=" * 78)

    exact = []

    for name, rows in CHANNELS.items():

        print()
        print(name)

        for order in (1, 2, 3):

            for degree in (0, 1, 2, 3):

                result = solve_operator(
                    rows,
                    order,
                    degree,
                )

                is_exact = (
                    result["unique"]
                    and verify_operator(
                        rows,
                        order,
                        degree,
                        result["solution"],
                    )
                )

                tag = (
                    "UNIQUE-EXACT"
                    if is_exact
                    else ""
                )

                print(
                    f"  order={order} "
                    f"basis_degree={degree}: "
                    f"consistent={result['consistent']} "
                    f"rank={result['rank']} "
                    f"unknowns={result['unknowns']} "
                    f"nullity={result['nullity']} "
                    f"{tag}"
                )

                if is_exact:
                    exact.append(
                        (
                            name,
                            order,
                            degree,
                            result["solution"],
                        )
                    )

    return exact


# -----------------------------------------------------------------------------
# LOCAL TRANSITION SEARCH
# -----------------------------------------------------------------------------

def local_search():

    print()
    print("=" * 78)
    print("3. LOCAL FALLING-BASIS TRANSITION SEARCH")
    print("=" * 78)

    for name, rows in CHANNELS.items():

        print()
        print(name)

        ps = sorted(rows)

        for p0, p1 in zip(ps, ps[1:]):

            subrows = {
                p0: rows[p0],
                p1: rows[p1],
            }

            print(
                f"  transition {p0}->{p1}"
            )

            for order in (1, 2, 3):

                for degree in (0, 1, 2):

                    result = solve_operator(
                        subrows,
                        order,
                        degree,
                    )

                    is_exact = (
                        result["unique"]
                        and verify_operator(
                            subrows,
                            order,
                            degree,
                            result["solution"],
                        )
                    )

                    if is_exact:
                        print(
                            f"    order={order} "
                            f"degree={degree} "
                            f"UNIQUE-EXACT"
                        )
                        print_operator(
                            result["solution"],
                            order,
                            degree,
                        )


# -----------------------------------------------------------------------------
# SECOND TEST:
# DOES A SIMPLE OPERATOR APPEAR AFTER RESCALING q[p,r] ?
# -----------------------------------------------------------------------------

def normalized_rows(rows):
    """
    Remove the first nonzero entry of each row.

    This is only a diagnostic normalization; exactness is retained.
    """
    out = {}

    for p, row in rows.items():

        nz = next(
            (x for x in row if x != 0),
            Fraction(1),
        )

        out[p] = [
            Fraction(x) / nz
            for x in row
        ]

    return out


def normalized_search():

    print()
    print("=" * 78)
    print("4. ROW-SCALE NORMALIZED FALLING-BASIS SEARCH")
    print("=" * 78)

    for name, rows in CHANNELS.items():

        norm = normalized_rows(rows)

        hits = []

        for order in (1, 2, 3):

            for degree in (0, 1, 2):

                result = solve_operator(
                    norm,
                    order,
                    degree,
                )

                exact = (
                    result["unique"]
                    and verify_operator(
                        norm,
                        order,
                        degree,
                        result["solution"],
                    )
                )

                if exact:
                    hits.append(
                        (order, degree)
                    )

        print(
            f"  {name}: "
            f"unique_exact_candidates={hits}"
        )


# -----------------------------------------------------------------------------
# RECONSTRUCTION
# -----------------------------------------------------------------------------

def reconstruction_sanity():

    print()
    print("=" * 78)
    print("5. EXACT RECONSTRUCTION SANITY")
    print("=" * 78)

    ok = True

    for name, rows in CHANNELS.items():

        print(name)

        for p, row in rows.items():

            D = len(row) - 1

            # Evaluate the falling polynomial on n=0,...,D.
            values = []

            for n in range(D + 1):

                value = sum(
                    coeff * falling(n, r)
                    for r, coeff in enumerate(row)
                )

                values.append(value)

            recovered = falling_basis_interpolate(values)

            exact = recovered == row
            ok &= exact

            print(
                f"  p={p}: exact={exact}"
            )

        print()

    return ok


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main():

    failures = 0

    print("=" * 78)
    print("EXPERIMENT 119 — EXACT FALLING-BASIS r-OPERATOR AUDIT")
    print("=" * 78)

    compare_bodd()
    exact_candidates = universal_search()
    local_search()
    normalized_search()

    reconstruction_ok = reconstruction_sanity()

    if not reconstruction_ok:
        failures += 1

    print()
    print("=" * 78)
    print("6. STRUCTURAL INTERPRETATION")
    print("=" * 78)
    print()
    print(
        "The monomial-r operator from Experiment 118 is exact, "
        "but its coefficients are extremely large."
    )
    print()
    print(
        "Experiment 119 changes only the basis:"
    )
    print()
    print(
        "    F_a(r) = sum_b c[a,b] r_(b)."
    )
    print()
    print(
        "The purpose is to determine whether the apparently "
        "complicated operator becomes sparse or small in the "
        "falling-factorial basis."
    )
    print()
    print(
        "A UNIQUE-EXACT operator is structurally meaningful."
    )
    print(
        "Underdetermined fits are deliberately not promoted."
    )
    print()
    print(
        "A particularly strong positive result would be:"
    )
    print(
        "    * lower operator order;"
    )
    print(
        "    * lower falling degree;"
    )
    print(
        "    * small integer/rational coefficients;"
    )
    print(
        "    * the same operator across multiple sectors."
    )
    print()
    print(
        "All computations use exact Fraction arithmetic."
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

    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)
    print(
        f"  reconstruction={reconstruction_ok}"
    )
    print(
        f"  universal_candidates={len(exact_candidates)}"
    )
    print(
        f"  failures={failures}"
    )
    print(
        f"  ALL BASIC CHECKS PASS={failures == 0}"
    )

    print()
    print("EXPERIMENT 119 COMPLETE")


if __name__ == "__main__":
    main()

