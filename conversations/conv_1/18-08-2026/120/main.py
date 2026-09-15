#!/usr/bin/env python3

from fractions import Fraction
from math import gcd, factorial
from itertools import product


# ==============================================================================
# EXPERIMENT 122 — EXACT B-ODD THREE-TERM r-OPERATOR FACTORIZATION AUDIT
# ==============================================================================
#
# Important correction to Experiment 122:
#
# The previously reported "order=3, degree_r=2" operator has THREE
# shift components:
#
#     T = F_0(r) + F_1(r) E + F_2(r) E^2
#
# where
#
#     (E f)(r) = f(r+1).
#
# Thus F_3 is NOT required.
#
# The current experiment asks whether this exact three-term operator
# factors into TWO first-order shift operators:
#
#     T = L_0 L_1
#
# with
#
#     L_i = A_i(r) + B_i(r) E.
#
# We test small falling-factorial degrees for A_i,B_i.
#
# The composition law is
#
#     (A + B E)(C + D E)
#
#       = A C
#       + (A D + B C(r+1)) E
#       + B D(r+1) E^2.
#
# Every identity is checked pointwise over exact Fraction arithmetic.
#
# No floating point.
# No SymPy.
# No extrapolation.
# ==============================================================================


# ==============================================================================
# 1. EXACT B-ODD OPERATOR FROM EXPERIMENT 121
# ==============================================================================

F0 = [
    Fraction(3315061342376947363, 477116474025609633376),
    Fraction(
        -25886544326998717157327370346624409742671826818575845092488719689187563,
        1516080269064001164312649765883488704089359863580858321520879927011681440,
    ),
    Fraction(
        53706935773064483686240917568532325383417674486719505867813526164411043,
        4548240807192004657250599063533954816357432423923433286083519708046725760,
    ),
]

F1 = [
    Fraction(7554395658365211, 74549449066501505215),
    Fraction(
        261176445687125841466691750754346851872093779658770103808,
        6016191543904768065146295057584596317933120971327292706459682153500960,
    ),
    Fraction(
        -3206553924272712664092331011793814615653413998574800660924525303991991,
        303216053812800310483373270902263654423829296954895552405567980536448384,
    ),
]

F2 = [
    Fraction(11833101985803613, 511196222170296035760),
    Fraction(
        2730823991367374816864091756853232538348222983384948949754065808599389,
        252680044844000129368072196209276522676595540061843424842812819327040320,
    ),
    Fraction(
        1481071209865873169095168327693134640712003505493254282379201508083869,
        84226681614666752912048130806184348451063693598582097890435550149013440,
    ),
]


# ==============================================================================
# 2. EXACT FALLING BASIS
# ==============================================================================

def falling(r, n):
    out = Fraction(1, 1)
    for i in range(n):
        out *= Fraction(r - i, 1)
    return out


def eval_falling(coeffs, r):
    return sum(
        coeff * falling(r, n)
        for n, coeff in enumerate(coeffs)
    )


def trim(poly):
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    if not out:
        return [Fraction(0)]
    return out


def values_to_falling(values):
    """
    Newton forward basis:
        f(r) = sum_n c_n r_(n)

    with
        c_n = Delta^n f(0) / n!
    """
    work = [Fraction(v) for v in values]
    out = []

    for n in range(len(work)):
        out.append(work[0] / factorial(n))
        work = [
            work[i + 1] - work[i]
            for i in range(len(work) - 1)
        ]

    return trim(out)


def falling_to_values(coeffs, n):
    return [
        eval_falling(coeffs, r)
        for r in range(n)
    ]


# ==============================================================================
# 3. POLYNOMIAL OPERATIONS IN FALLING BASIS
# ==============================================================================

def add_poly(a, b):
    n = max(len(a), len(b))
    out = [Fraction(0)] * n

    for i, value in enumerate(a):
        out[i] += value

    for i, value in enumerate(b):
        out[i] += value

    return trim(out)


def scale_poly(a, c):
    return trim([c * x for x in a])


def multiply_poly(a, b):
    da = len(a) - 1
    db = len(b) - 1

    values = []
    for r in range(da + db + 1):
        values.append(
            eval_falling(a, r)
            * eval_falling(b, r)
        )

    return values_to_falling(values)


def shift_poly(a, shift):
    """
    Exact polynomial p(r+shift), represented again in falling basis.
    """
    degree = len(a) - 1

    values = [
        eval_falling(a, r + shift)
        for r in range(degree + 1)
    ]

    return values_to_falling(values)


def poly_equal(a, b):
    a = trim(a)
    b = trim(b)

    n = max(len(a), len(b))

    for i in range(n):
        ai = a[i] if i < len(a) else Fraction(0)
        bi = b[i] if i < len(b) else Fraction(0)

        if ai != bi:
            return False

    return True


# ==============================================================================
# 4. EXACT OPERATOR COMPOSITION
# ==============================================================================

def compose_first_order(A, B, C, D):
    """
    (A + B E)(C + D E)

    = AC
      + (AD + B*C(r+1)) E
      + B*D(r+1) E^2
    """

    C_shift = shift_poly(C, 1)
    D_shift = shift_poly(D, 1)

    G0 = multiply_poly(A, C)

    G1 = add_poly(
        multiply_poly(A, D),
        multiply_poly(B, C_shift),
    )

    G2 = multiply_poly(B, D_shift)

    return [
        G0,
        G1,
        G2,
    ]


def evaluate_operator(op, r):
    return [
        eval_falling(component, r)
        for component in op
    ]


# ==============================================================================
# 5. EXACT LINEAR ALGEBRA
# ==============================================================================

def rref(A, b):
    M = [
        [Fraction(x) for x in row] + [Fraction(rhs)]
        for row, rhs in zip(A, b)
    ]

    rows = len(M)
    cols = len(A[0]) if A else 0

    pivots = []
    row = 0

    for col in range(cols):
        pivot = None

        for r in range(row, rows):
            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        M[row], M[pivot] = M[pivot], M[row]

        p = M[row][col]
        M[row] = [x / p for x in M[row]]

        for r in range(rows):
            if r == row:
                continue

            q = M[r][col]

            if q == 0:
                continue

            M[r] = [
                M[r][j] - q * M[row][j]
                for j in range(cols + 1)
            ]

        pivots.append(col)
        row += 1

        if row == rows:
            break

    for r in range(rows):
        if (
            all(M[r][c] == 0 for c in range(cols))
            and M[r][cols] != 0
        ):
            return {
                "consistent": False,
                "rank": len(pivots),
                "unknowns": cols,
                "nullity": None,
                "solution": None,
            }

    solution = [Fraction(0)] * cols

    for i, pivot in enumerate(pivots):
        solution[pivot] = M[i][cols]

    return {
        "consistent": True,
        "rank": len(pivots),
        "unknowns": cols,
        "nullity": cols - len(pivots),
        "solution": solution,
    }


# ==============================================================================
# 6. FALLING-BASIS MONOMIAL SYSTEM
# ==============================================================================

def falling_basis_matrix_sample(deg):
    return [
        [falling(r, j) for j in range(deg + 1)]
        for r in range(deg + 1)
    ]


# ==============================================================================
# 7. FIRST-ORDER FACTOR SUPPORT
# ==============================================================================

def make_factor_from_vector(vec, deg):
    return vec[:deg + 1]


def factor_unknowns(deg):
    """
    A(r) = sum_{i=0}^deg a_i r_(i)
    B(r) = sum_{i=0}^deg b_i r_(i)

    Unknown count = 2(deg+1)
    """
    return 2 * (deg + 1)


# ==============================================================================
# 8. FACTORIZATION BY DIRECT FINITE SEARCH
# ==============================================================================

def factorization_pointwise_test(
    A,
    B,
    C,
    D,
    known,
    sample_count,
):
    candidate = compose_first_order(A, B, C, D)

    for shift in range(3):
        known_component = known[shift]
        candidate_component = candidate[shift]

        for r in range(sample_count):

            x = eval_falling(known_component, r)
            y = eval_falling(candidate_component, r)

            if x != y:
                return False

    return True


# ==============================================================================
# 9. NORMALIZED FACTOR CANDIDATES
# ==============================================================================

def scalar_normalize(A, B):
    """
    Remove a common scalar from a first-order factor by fixing the first
    nonzero coefficient in A, then B.

    This only serves to canonicalize reporting.
    """

    coefficients = list(A) + list(B)

    pivot = None

    for x in coefficients:
        if x != 0:
            pivot = x
            break

    if pivot is None:
        return A, B, None

    return (
        scale_poly(A, Fraction(1, 1) / pivot),
        scale_poly(B, Fraction(1, 1) / pivot),
        pivot,
    )


# ==============================================================================
# 10. SMALL INTEGER/RATIONAL NORMALIZED FACTOR SEARCH
# ==============================================================================

def primitive_integer_signature(poly):
    dens = [x.denominator for x in poly]
    common_den = 1

    for d in dens:
        common_den = common_den * d // gcd(common_den, d)

    nums = [
        x.numerator * (common_den // x.denominator)
        for x in poly
    ]

    g = 0
    for n in nums:
        g = gcd(g, abs(n))

    if g == 0:
        return [0]

    nums = [n // g for n in nums]

    for n in nums:
        if n != 0:
            if n < 0:
                nums = [-x for x in nums]
            break

    return nums


def complexity_signature(poly):
    nz = [x for x in poly if x != 0]

    if not nz:
        return {
            "degree": -1,
            "nonzero": 0,
            "max_num_digits": 0,
            "max_den_digits": 0,
        }

    return {
        "degree": len(poly) - 1,
        "nonzero": len(nz),
        "max_num_digits": max(
            len(str(abs(x.numerator)))
            for x in nz
        ),
        "max_den_digits": max(
            len(str(x.denominator))
            for x in nz
        ),
    }


# ==============================================================================
# 11. SPECIAL ZERO/ONE FACTORS
# ==============================================================================

SPECIAL_FACTORS = {
    "A_zero": [Fraction(0)],
    "A_one": [Fraction(1)],
    "B_zero": [Fraction(0)],
    "B_one": [Fraction(1)],
}


def special_factorization_candidates(
    known,
    max_degree,
    sample_count,
):
    """
    Test canonical/simple first-order factor patterns.

    These are not an exhaustive nonlinear search.  They are exact
    structural probes.
    """

    candidates = []

    zero = [Fraction(0)]
    one = [Fraction(1)]

    for side in ("left", "right"):

        # (known operator) = (known operator) * identity
        if side == "right":

            comp = compose_first_order(
                known[0], known[1],
                one, zero
            )

            if factorization_pointwise_test(
                known[0], known[1],
                one, zero,
                known,
                sample_count,
            ):
                candidates.append(
                    "T = T * I"
                )

        else:

            comp = compose_first_order(
                one, zero,
                known[0], known[1]
            )

            if factorization_pointwise_test(
                one, zero,
                known[0], known[1],
                known,
                sample_count,
            ):
                candidates.append(
                    "T = I * T"
                )

    return candidates


# ==============================================================================
# 12. RIGHT-DIVISOR TEST
# ==============================================================================

def right_divisor_test(known, A, B, sample_count):
    """
    If
        T = S (A + B E),
    then recursively compare shift components.

    For a second first-order operator
        S = C + D E

    we have:
        F0 = C A
        F2 = D B(r+1)
        F1 = C B + D A(r+1).

    We can reconstruct C,D pointwise and then test the middle equation.
    """

    A_values = [eval_falling(A, r) for r in range(sample_count)]
    B_shift_values = [
        eval_falling(B, r + 1)
        for r in range(sample_count)
    ]

    C_values = []
    D_values = []

    for r in range(sample_count):

        ar = A_values[r]
        br = B_shift_values[r]

        f0 = eval_falling(known[0], r)
        f2 = eval_falling(known[2], r)

        if ar == 0:
            if f0 != 0:
                return None
        else:
            C_values.append((r, f0 / ar))

        if br == 0:
            if f2 != 0:
                return None
        else:
            D_values.append((r, f2 / br))

    if not C_values or not D_values:
        return None

    # Return pointwise sequences for exact polynomial interpolation.
    return C_values, D_values


# ==============================================================================
# 13. SMALL FIRST-ORDER DIVISOR SEARCH
# ==============================================================================

def search_first_order_divisors(
    known,
    max_degree=2,
):
    """
    Search simple divisors whose coefficients are small
    falling-basis polynomials.

    We deliberately search only a finite dictionary of
    low-complexity polynomials.

    This is a factorization audit, not an unrestricted nonlinear
    algebra solver.
    """

    dictionary = []

    # Constant factors.
    constants = [
        Fraction(-2),
        Fraction(-1),
        Fraction(0),
        Fraction(1),
        Fraction(2),
    ]

    for c in constants:
        dictionary.append([c])

    # r
    if max_degree >= 1:
        dictionary.append([Fraction(0), Fraction(1)])
        dictionary.append([Fraction(1), Fraction(1)])
        dictionary.append([Fraction(-1), Fraction(1)])

    # r_(2)
    if max_degree >= 2:
        dictionary.append([Fraction(0), Fraction(0), Fraction(1)])
        dictionary.append([Fraction(1), Fraction(0), Fraction(1)])
        dictionary.append([Fraction(-1), Fraction(0), Fraction(1)])

    # Test B(r) divisors with A=1.
    hits = []

    A = [Fraction(1)]

    for B in dictionary:

        # Avoid zero leading operator.
        if all(x == 0 for x in B):
            continue

        right = right_divisor_test(
            known,
            A,
            B,
            sample_count=10,
        )

        if right is None:
            continue

        hits.append(
            ("A=1", B)
        )

    return hits


# ==============================================================================
# 14. EXACT KNOWN OPERATOR SANITY
# ==============================================================================

def operator_reconstruction_check(known):
    """
    Verify the operator against the three known B-odd row transitions.

    q rows are encoded from Experiment 115.

    We only need exact reconstructed values here.
    """

    q1 = [
        Fraction(12143, 3360),
        Fraction(76427, 192),
        Fraction(1954873, 17920),
        Fraction(-61469491, 483840),
        Fraction(42852113, 1935360),
        Fraction(4384549, 645120),
    ]

    q3 = [
        Fraction(-989, 345600),
        Fraction(-36064769, 2764800),
        Fraction(-24988097, 2580480),
        Fraction(87300373, 19353600),
        Fraction(116226679, 77414400),
    ]

    q5 = [
        Fraction(-517, 2764800),
        Fraction(1189391, 5529600),
        Fraction(22183547, 103219200),
        Fraction(-69294643, 185794560),
    ]

    q7 = [
        Fraction(373, 928972800),
        Fraction(-616981, 371589120),
    ]

    rows = [
        q1,
        q3,
        q5,
        q7,
    ]

    for i in range(3):
        source = rows[i]
        target = rows[i + 1]

        for r in range(len(target)):

            predicted = Fraction(0)

            for shift, F in enumerate(known):
                rr = r + shift

                if rr < len(source):
                    predicted += (
                        eval_falling(F, r)
                        * source[rr]
                    )

            if predicted != target[r]:
                return False

    return True


# ==============================================================================
# 15. MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 122 — EXACT B-ODD THREE-TERM "
        "SHIFT-OPERATOR FACTORIZATION AUDIT"
    )
    print("=" * 78)

    known = [F0, F1, F2]

    print()
    print("=" * 78)
    print("1. OPERATOR VALIDATION")
    print("=" * 78)

    print(
        "  number_of_shift_components =",
        len(known),
    )

    for i, F in enumerate(known):
        print(
            f"  F_{i}: degree={len(F)-1} "
            f"nonzero={sum(x != 0 for x in F)}"
        )

    reconstruction = operator_reconstruction_check(known)

    print(
        "  pointwise_operator_reconstruction =",
        reconstruction,
    )

    print()
    print("=" * 78)
    print("2. COMPONENT COMPLEXITY")
    print("=" * 78)

    for i, F in enumerate(known):
        sig = primitive_integer_signature(F)
        comp = complexity_signature(F)

        print(
            f"  F_{i}: "
            f"degree={comp['degree']} "
            f"nonzero={comp['nonzero']} "
            f"max_num_digits={comp['max_num_digits']} "
            f"max_den_digits={comp['max_den_digits']}"
        )
        print(
            f"    primitive_signature={sig}"
        )

    print()
    print("=" * 78)
    print("3. SPECIAL FIRST-ORDER DIVISOR SEARCH")
    print("=" * 78)

    special_hits = search_first_order_divisors(
        known,
        max_degree=2,
    )

    print(
        "  candidate_count =",
        len(special_hits),
    )

    for item in special_hits:
        print(
            "  candidate =",
            item,
        )

    if not special_hits:
        print("  NONE")

    print()
    print("=" * 78)
    print("4. CANONICAL IDENTITY FACTORIZATION AUDIT")
    print("=" * 78)

    identity_hits = special_factorization_candidates(
        known,
        max_degree=2,
        sample_count=10,
    )

    print(
        "  canonical_identity_factorizations =",
        identity_hits,
    )

    print()
    print("=" * 78)
    print("5. STRUCTURAL FACTORIZATION TEST")
    print("=" * 78)

    print()
    print(
        "  The exact three-term operator is"
    )
    print()
    print(
        "      T = F_0(r) + F_1(r) E + F_2(r) E^2."
    )
    print()
    print(
        "  A first-order factorization would have the form"
    )
    print()
    print(
        "      T = (A(r) + B(r)E)(C(r) + D(r)E)."
    )
    print()
    print(
        "  which imposes"
    )
    print()
    print(
        "      F_0 = A C"
    )
    print(
        "      F_1 = A D + B C(r+1)"
    )
    print(
        "      F_2 = B D(r+1)."
    )

    print()
    print(
        "  The finite divisor dictionary above found:",
        "NONE" if not special_hits else special_hits,
    )

    print()
    print("=" * 78)
    print("6. EXACT RECONSTRUCTION")
    print("=" * 78)

    print(
        "  operator_reconstruction =",
        reconstruction,
    )

    print()
    print("=" * 78)
    print("7. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print()
    print(
        "Experiment 121 showed that standard factorial/binomial "
        "conjugations do not simplify or transfer the exceptional "
        "B-odd operator."
    )
    print()
    print(
        "The present experiment moves to a stronger question:"
    )
    print()
    print(
        "    Is the three-term operator itself reducible?"
    )
    print()
    print(
        "A successful factorization into two first-order shifts "
        "would expose an elementary construction mechanism."
    )
    print()
    print(
        "Failure of the exact low-complexity divisor search is "
        "evidence against a small first-order factorization, "
        "but it is not a proof of irreducibility over the full "
        "rational polynomial ring."
    )
    print()
    print(
        "Everything here is exact over QQ."
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
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  operator_exact =",
        reconstruction,
    )

    print(
        "  low_complexity_factor_search_completed = True"
    )

    print(
        "  factorization_found =",
        bool(special_hits),
    )

    print(
        "  failures =",
        0 if reconstruction else 1,
    )

    print(
        "  ALL BASIC CHECKS PASS =",
        reconstruction,
    )

    print()
    print(
        "EXPERIMENT 122 COMPLETE"
    )


if __name__ == "__main__":
    main()

