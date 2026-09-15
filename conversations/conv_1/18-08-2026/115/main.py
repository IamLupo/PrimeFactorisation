#!/usr/bin/env python3

from fractions import Fraction
from math import gcd


# =============================================================================
# EXPERIMENT 118 — EXACT UNIVERSAL p/r OPERATOR + FACTORIZATION AUDIT
# FIXED VERSION
# =============================================================================
#
# Fixes:
#   * imports gcd explicitly from math;
#   * keeps all arithmetic exact with Fraction;
#   * avoids SymPy entirely;
#   * handles zero coefficient blocks safely;
#   * preserves the Experiment 117 operator search.
#
# Main recurrence:
#
#     q[p+2,r] = sum_{a=0}^{order-1} F_a(p,r) q[p,r+a]
#
# where
#
#     F_a(p,r)
#       = sum_{i=0}^dp sum_{j=0}^dr c[a,i,j] p^i r^j.
#
# A solution is promoted only when it is:
#   1. consistent,
#   2. unique,
#   3. verified pointwise on every observed transition.
#
# No floating point.
# No SymPy.
# No extrapolation.
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
    """Exact reduced row echelon form over Fraction."""
    if not matrix:
        return [], [], 0

    M = [
        [Fraction(x) for x in row]
        for row in matrix
    ]

    rows = len(M)
    cols = len(M[0])

    pivots = []
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

        pv = M[pivot_row][col]
        M[pivot_row] = [
            x / pv for x in M[pivot_row]
        ]

        for r in range(rows):
            if r == pivot_row:
                continue

            factor = M[r][col]

            if factor == 0:
                continue

            M[r] = [
                M[r][j] - factor * M[pivot_row][j]
                for j in range(cols)
            ]

        pivots.append(col)
        pivot_row += 1

        if pivot_row == rows:
            break

    return M, pivots, pivot_row


def rank(matrix):
    if not matrix:
        return 0
    return rref(matrix)[2]


def solve_unique(A, b):
    """
    Solve A x = b exactly.

    Returns:
        consistent
        unique
        rank
        nullity
        solution
    """
    if not A:
        return {
            "consistent": False,
            "unique": False,
            "rank": 0,
            "nullity": None,
            "solution": None,
        }

    nvars = len(A[0])

    augmented = [
        [Fraction(x) for x in row] + [Fraction(bi)]
        for row, bi in zip(A, b)
    ]

    R, pivots, rr = rref(augmented)

    matrix_rank = rank(A)

    # Inconsistent row.
    for row in R:
        if all(row[j] == 0 for j in range(nvars)) and row[nvars] != 0:
            return {
                "consistent": False,
                "unique": False,
                "rank": matrix_rank,
                "nullity": None,
                "solution": None,
            }

    nullity = nvars - matrix_rank

    if nullity != 0:
        return {
            "consistent": True,
            "unique": False,
            "rank": matrix_rank,
            "nullity": nullity,
            "solution": None,
        }

    solution = [Fraction(0)] * nvars

    for i, pivot in enumerate(pivots):
        if pivot < nvars:
            solution[pivot] = R[i][nvars]

    return {
        "consistent": True,
        "unique": True,
        "rank": matrix_rank,
        "nullity": 0,
        "solution": solution,
    }


# -----------------------------------------------------------------------------
# FALLING FACTORIAL
# -----------------------------------------------------------------------------

def falling(n, r):
    if r < 0 or r > n:
        return Fraction(0)

    out = 1

    for i in range(r):
        out *= n - i

    return Fraction(out)


# -----------------------------------------------------------------------------
# p/r MONOMIAL BASIS
# -----------------------------------------------------------------------------

def monomial_basis(p, r, degree_p, degree_r):
    values = []

    for ip in range(degree_p + 1):
        p_power = p ** ip

        for ir in range(degree_r + 1):
            values.append(
                Fraction(p_power * (r ** ir))
            )

    return values


# -----------------------------------------------------------------------------
# BUILD UNIVERSAL p/r OPERATOR SYSTEM
# -----------------------------------------------------------------------------

def build_pr_system(rows, order, degree_p, degree_r):
    """
    Build equations for

        q[p_next,r] =
            sum_a F_a(p,r) q[p,r+a].

    Each F_a has bidegree
        degree_p in p
        degree_r in r.
    """
    ps = sorted(rows)

    block_size = (
        (degree_p + 1)
        * (degree_r + 1)
    )

    total_unknowns = order * block_size
    equations = []

    for p0, p1 in zip(ps, ps[1:]):
        prev = rows[p0]
        nxt = rows[p1]

        # Only use the actual target support.
        for r in range(len(nxt)):

            lhs_row = [
                Fraction(0)
            ] * total_unknowns

            for a in range(order):
                source_index = r + a

                if source_index >= len(prev):
                    continue

                source = prev[source_index]

                basis = monomial_basis(
                    p0,
                    r,
                    degree_p,
                    degree_r,
                )

                offset = a * block_size

                for j, basis_value in enumerate(basis):
                    lhs_row[offset + j] += (
                        source * basis_value
                    )

            # Operator expression - target = 0.
            equations.append(
                lhs_row + [-nxt[r]]
            )

    return equations


def solve_pr_operator(
    rows,
    order,
    degree_p,
    degree_r,
):
    system = build_pr_system(
        rows,
        order,
        degree_p,
        degree_r,
    )

    if not system:
        unknowns = (
            order
            * (degree_p + 1)
            * (degree_r + 1)
        )

        return {
            "consistent": False,
            "unique": False,
            "rank": 0,
            "unknowns": unknowns,
            "nullity": None,
            "solution": None,
        }

    A = [
        row[:-1]
        for row in system
    ]

    b = [
        -row[-1]
        for row in system
    ]

    result = solve_unique(A, b)
    result["unknowns"] = len(A[0])

    return result


# -----------------------------------------------------------------------------
# OPERATOR EVALUATION
# -----------------------------------------------------------------------------

def eval_operator(
    prev,
    solution,
    p,
    r,
    order,
    degree_p,
    degree_r,
):
    block_size = (
        (degree_p + 1)
        * (degree_r + 1)
    )

    total = Fraction(0)

    for a in range(order):
        source_index = r + a

        if source_index >= len(prev):
            continue

        q = prev[source_index]

        block_start = a * block_size

        coeff = Fraction(0)
        offset = 0

        for ip in range(degree_p + 1):
            p_power = p ** ip

            for ir in range(degree_r + 1):
                coeff += (
                    solution[block_start + offset]
                    * p_power
                    * (r ** ir)
                )
                offset += 1

        total += coeff * q

    return total


def verify_operator(
    rows,
    order,
    degree_p,
    degree_r,
    solution,
):
    ps = sorted(rows)

    for p0, p1 in zip(ps, ps[1:]):
        prev = rows[p0]
        nxt = rows[p1]

        for r in range(len(nxt)):
            predicted = eval_operator(
                prev,
                solution,
                p0,
                r,
                order,
                degree_p,
                degree_r,
            )

            if predicted != nxt[r]:
                return False

    return True


# -----------------------------------------------------------------------------
# COEFFICIENT-BLOCK UTILITIES
# -----------------------------------------------------------------------------

def split_blocks(solution, order, degree_p, degree_r):
    block_size = (
        (degree_p + 1)
        * (degree_r + 1)
    )

    return [
        solution[
            a * block_size:
            (a + 1) * block_size
        ]
        for a in range(order)
    ]


def common_integer_factor(coeffs):
    """
    Return gcd of numerators divided by lcm of denominators.

    FIX:
      * imports math.gcd;
      * handles a single nonzero coefficient;
      * handles all-zero coefficient vectors.
    """
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

        denominator_lcm = (
            denominator_lcm
            // gcd(denominator_lcm, c.denominator)
            * c.denominator
        )

    return Fraction(
        numerator_gcd,
        denominator_lcm,
    )


def polynomial_terms(coeffs, degree_p, degree_r):
    terms = []
    idx = 0

    for ip in range(degree_p + 1):
        for ir in range(degree_r + 1):

            c = coeffs[idx]
            idx += 1

            if c == 0:
                continue

            factors = []

            if ip == 1:
                factors.append("p")
            elif ip > 1:
                factors.append(f"p^{ip}")

            if ir == 1:
                factors.append("r")
            elif ir > 1:
                factors.append(f"r^{ir}")

            monomial = "*".join(factors)

            if monomial:
                if c == 1:
                    term = monomial
                elif c == -1:
                    term = "-" + monomial
                else:
                    term = f"({c})*{monomial}"
            else:
                term = str(c)

            terms.append(term)

    return terms


def polynomial_string(coeffs, degree_p, degree_r):
    terms = polynomial_terms(
        coeffs,
        degree_p,
        degree_r,
    )

    if not terms:
        return "0"

    return " + ".join(terms)


def print_operator(
    solution,
    order,
    degree_p,
    degree_r,
):
    blocks = split_blocks(
        solution,
        order,
        degree_p,
        degree_r,
    )

    for a, coeffs in enumerate(blocks):
        print(
            f"    F_{a}(p,r) = "
            f"{polynomial_string(coeffs, degree_p, degree_r)}"
        )


# -----------------------------------------------------------------------------
# SPECIAL B-ODD AUDIT
# -----------------------------------------------------------------------------

def audit_bodd_special():
    print("=" * 78)
    print("SPECIAL B-ODD OPERATOR AUDIT")
    print("=" * 78)

    rows = B_odd

    order = 3
    degree_p = 0
    degree_r = 2

    result = solve_pr_operator(
        rows,
        order,
        degree_p,
        degree_r,
    )

    print(
        f"  order={order} "
        f"degree_p={degree_p} "
        f"degree_r={degree_r}"
    )
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

    if not result["unique"]:
        print("  operator not uniquely determined")
        return False

    exact = verify_operator(
        rows,
        order,
        degree_p,
        degree_r,
        result["solution"],
    )

    print(
        f"  pointwise_exact={exact}"
    )

    if exact:
        print_operator(
            result["solution"],
            order,
            degree_p,
            degree_r,
        )

        print("  coefficient-block common factors:")

        blocks = split_blocks(
            result["solution"],
            order,
            degree_p,
            degree_r,
        )

        for a, coeffs in enumerate(blocks):
            print(
                f"    F_{a}: "
                f"{common_integer_factor(coeffs)}"
            )

    return exact


# -----------------------------------------------------------------------------
# SYSTEMATIC SEARCH
# -----------------------------------------------------------------------------

def systematic_search():
    print()
    print("=" * 78)
    print("SYSTEMATIC UNIVERSAL p/r SEARCH")
    print("=" * 78)

    exact_candidates = []

    for name, rows in CHANNELS.items():

        print()
        print(name)

        for order in (1, 2, 3):

            for degree_p in (0, 1, 2):

                for degree_r in (0, 1, 2):

                    result = solve_pr_operator(
                        rows,
                        order,
                        degree_p,
                        degree_r,
                    )

                    exact = False

                    if result["unique"]:
                        exact = verify_operator(
                            rows,
                            order,
                            degree_p,
                            degree_r,
                            result["solution"],
                        )

                    tag = (
                        "UNIQUE-EXACT"
                        if exact
                        else ""
                    )

                    print(
                        f"  order={order} "
                        f"dp={degree_p} "
                        f"dr={degree_r}: "
                        f"consistent={result['consistent']} "
                        f"rank={result['rank']} "
                        f"unknowns={result['unknowns']} "
                        f"nullity={result['nullity']} "
                        f"{tag}"
                    )

                    if exact:
                        exact_candidates.append(
                            (
                                name,
                                order,
                                degree_p,
                                degree_r,
                                result["solution"],
                            )
                        )

    return exact_candidates


# -----------------------------------------------------------------------------
# LOWEST-COMPLEXITY UNIQUE OPERATORS
# -----------------------------------------------------------------------------

def print_lowest_candidates(candidates):
    print()
    print("=" * 78)
    print("LOWEST-COMPLEXITY UNIQUE EXACT OPERATORS")
    print("=" * 78)

    for name in CHANNELS:

        channel_candidates = [
            item
            for item in candidates
            if item[0] == name
        ]

        if not channel_candidates:
            print(f"{name}: NONE")
            continue

        channel_candidates.sort(
            key=lambda item: (
                item[1],                    # order
                item[2] + item[3],          # total degree
                item[2],
                item[3],
            )
        )

        (
            cname,
            order,
            degree_p,
            degree_r,
            solution,
        ) = channel_candidates[0]

        print()
        print(
            f"{cname}: "
            f"order={order}, "
            f"degree_p={degree_p}, "
            f"degree_r={degree_r}"
        )

        print_operator(
            solution,
            order,
            degree_p,
            degree_r,
        )


# -----------------------------------------------------------------------------
# RECONSTRUCTION SANITY
# -----------------------------------------------------------------------------

def reconstruction_sanity():
    print()
    print("=" * 78)
    print("EXACT FALLING-BASIS RECONSTRUCTION SANITY")
    print("=" * 78)

    overall = True

    for name, rows in CHANNELS.items():

        print(name)

        for p, row in rows.items():

            D = len(row) - 1

            values = []

            for n in range(D + 1):
                total = Fraction(0)

                for r, coeff in enumerate(row):
                    total += (
                        coeff
                        * falling(n, r)
                    )

                values.append(total)

            recovered = [
                Fraction(0)
            ] * (D + 1)

            for n in range(D + 1):

                residual = values[n]

                for r in range(n):
                    residual -= (
                        recovered[r]
                        * falling(n, r)
                    )

                recovered[n] = (
                    residual
                    / falling(n, n)
                )

            exact = recovered == row
            overall &= exact

            print(
                f"  p={p}: exact={exact}"
            )

        print()

    return overall


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main():

    failures = 0

    print("=" * 78)
    print("EXPERIMENT 118 — EXACT UNIVERSAL p/r OPERATOR")
    print("AND FACTORIZATION AUDIT")
    print("=" * 78)
    print()

    # 1. Profile
    print("=" * 78)
    print("1. RESIDUAL ROW PROFILE")
    print("=" * 78)

    for name, rows in CHANNELS.items():
        print(name)

        for p in sorted(rows):
            print(
                f"  p={p}: degree={len(rows[p]) - 1}"
            )

        print()

    # 2. Special B-odd operator
    special_ok = audit_bodd_special()

    if not special_ok:
        failures += 1

    # 3. Systematic search
    candidates = systematic_search()

    # 4. Lowest candidates
    print_lowest_candidates(candidates)

    # 5. Reconstruction
    reconstruction_ok = reconstruction_sanity()

    if not reconstruction_ok:
        failures += 1

    # 6. Final interpretation
    print()
    print("=" * 78)
    print("6. STRUCTURAL INTERPRETATION")
    print("=" * 78)
    print()
    print(
        "Experiment 117 isolated one exact universal result:"
    )
    print()
    print(
        "  B-odd:"
    )
    print(
        "      order=3"
    )
    print(
        "      degree_p=0"
    )
    print(
        "      degree_r=2"
    )
    print()
    print(
        "This run independently reconstructs that operator"
    )
    print(
        "and checks whether comparable unique exact operators"
    )
    print(
        "exist in the other sectors at the same small complexity."
    )
    print()
    print(
        "Only unique exact operators are treated as structural."
    )
    print(
        "Underdetermined local fits are not promoted."
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

    # 7. Final
    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  B_odd_special_operator={special_ok}"
    )
    print(
        f"  reconstruction={reconstruction_ok}"
    )
    print(
        "  universal_search_completed=True"
    )
    print(
        f"  failures={failures}"
    )
    print(
        f"  ALL BASIC CHECKS PASS={failures == 0}"
    )

    print()
    print("EXPERIMENT 118 COMPLETE")


if __name__ == "__main__":
    main()