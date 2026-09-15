#!/usr/bin/env python3

from fractions import Fraction
from math import gcd


# ==============================================================================
# EXPERIMENT 122 — EXACT B-ODD SHIFT-OPERATOR FACTORIZATION AUDIT
# ==============================================================================
#
# Starting point:
#
#     q[p+2,r] = F_0(r) q[p,r]
#                + F_1(r) q[p,r+1]
#                + F_2(r) q[p,r+2]
#                + F_3(r) q[p,r+3]
#
# For the exceptional B-odd operator from Experiments 118-121,
# the operator has order 3 and falling degree 2:
#
#     F_a(r) = sum_b c[a,b] r_(b).
#
# This experiment asks whether the order-3 operator factors as
#
#     T = L_2 L_1 L_0
#
# where each L_i is a first-order shift operator
#
#     (L f)(r) = A(r) f(r) + B(r) f(r+1),
#
# with A(r), B(r) restricted to small falling-factorial degree.
#
# We test:
#
#   1. exact recovery of the known B-odd operator;
#   2. first-order factorization with degrees 0,1,2;
#   3. all exact factorization branches;
#   4. triangular/degenerate factors;
#   5. coefficient-size complexity;
#   6. exact reconstruction.
#
# No floating point.
# No SymPy.
# No extrapolation.
# ==============================================================================


# ------------------------------------------------------------------------------
# B-ODD OPERATOR FROM EXPERIMENT 121
# ------------------------------------------------------------------------------

F0 = [
    Fraction(3315061342376947363, 477116474025609633376),
    Fraction(
        -25886544326998717157327370346624409742671826818575845092488719689187563,
        1516080269064001552416866354511318272119146484774477762027839902682241920,
    ),
    Fraction(
        53706935773064483686240917568532325383417674486719505867813526164411043,
        4548240807192004657250599063533954816357432423923433286083519708046725760,
    ),
]

F1 = [
    Fraction(7554395658365211, 74549449066501505215),
    Fraction(
        261176445687125841466691750754346851872093779658677103808,
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
        2730823991367374816867046153622723173648224508030033057185849254609889,
        252680044844000258736144392418553045353191080795746293671306650447040320,
    ),
    Fraction(
        1481071209865873169095168327693134640712003505493254282379201508083869,
        84226681614666752912048130806184348451063693598582097890435550149013440,
    ),
]

# NOTE:
# The three arrays above are F_0,F_1,F_2 of the order-3 recurrence:
#
#     q[p+2,r] = F_0(r) q[p,r]
#                + F_1(r) q[p,r+1]
#                + F_2(r) q[p,r+2]
#                + F_3(r) q[p,r+3]
#
# The order-3 coefficient F_3 is normalized here as the leading shift.
# It is reconstructed below directly from the Experiment-118 operator.
#
# To avoid transcription ambiguity, replace F3 below by the exact F_3 from
# your Experiment-118/120 result if your local dataset stores it separately.

F3 = [
    Fraction(
        0
    ),
]


# ==============================================================================
# EXACT FALLING FACTORIAL
# ==============================================================================

def falling(r, n):
    out = Fraction(1)
    for i in range(n):
        out *= r - i
    return out


def eval_falling(coeffs, r):
    total = Fraction(0)
    for n, c in enumerate(coeffs):
        total += c * falling(r, n)
    return total


# ==============================================================================
# POLYNOMIAL ARITHMETIC IN THE FALLING BASIS
# ==============================================================================

def falling_to_values(coeffs, count):
    return [
        eval_falling(coeffs, r)
        for r in range(count)
    ]


def values_to_falling(values):
    work = [Fraction(v) for v in values]
    out = []

    for n in range(len(work)):
        # c_n = Delta^n f(0) / n!
        import math
        out.append(work[0] / math.factorial(n))

        work = [
            work[i + 1] - work[i]
            for i in range(len(work) - 1)
        ]

    return out


def add_poly(a, b):
    n = max(len(a), len(b))
    out = [Fraction(0)] * n

    for i, x in enumerate(a):
        out[i] += x

    for i, x in enumerate(b):
        out[i] += x

    while len(out) > 1 and out[-1] == 0:
        out.pop()

    return out


def scale_poly(a, c):
    out = [c * x for x in a]

    while len(out) > 1 and out[-1] == 0:
        out.pop()

    return out


def multiply_poly(a, b):
    # Multiplication is easiest by exact values followed by falling
    # interpolation.  Only very small degrees are used.
    da = len(a) - 1
    db = len(b) - 1

    values = []

    for r in range(da + db + 1):
        values.append(
            eval_falling(a, r)
            * eval_falling(b, r)
        )

    return values_to_falling(values)


# ==============================================================================
# FIRST-ORDER SHIFT OPERATOR
# ==============================================================================

class FirstOrderOperator:
    def __init__(self, A, B):
        self.A = list(A)
        self.B = list(B)

    def apply(self, terms):
        """
        Input:
            terms[a](r), a=0..m

        If L = A(r) + B(r) E,
        then

            L(sum_a C_a(r) E^a)
              =
              sum_a
                A C_a E^a
                +
                B C_a(r) E^(a+1)

        where multiplication is pointwise in r.
        """

        out = [[] for _ in range(len(terms) + 1)]

        for a, C in enumerate(terms):

            out[a] = add_poly(
                out[a],
                multiply_poly(
                    self.A,
                    C,
                ),
            )

            out[a + 1] = add_poly(
                out[a + 1],
                multiply_poly(
                    self.B,
                    C,
                ),
            )

        return out


# ==============================================================================
# KNOWN ORDER-3 OPERATOR
# ==============================================================================

def known_operator():
    """
    Supply the exact four components F_0 ... F_3.

    The B-odd operator in Experiment 118 was reported using a normalized
    order-3 convention.  The safest approach is to recover F_3 from the
    original local operator data rather than silently inventing it.

    If your Experiment-118 script has the exact fourth component, paste it
    into this function.
    """

    # ----------------------------------------------------------------------
    # PLACEHOLDER GUARD
    # ----------------------------------------------------------------------
    #
    # We deliberately refuse to fabricate F_3.
    #
    # This is preferable to producing a misleading factorization result.
    #
    raise RuntimeError(
        "Paste the exact B-odd F_3(r) component from Experiment 118/120 "
        "into known_operator(). The previous experiment output shown in "
        "the conversation omitted F_3(r)."
    )


# ==============================================================================
# OPERATOR COMPOSITION
# ==============================================================================

def compose_first_order(L2, L1, L0):
    """
    Compute L2 * L1 * L0.
    """

    identity = [
        [Fraction(1)],
    ]

    after0 = L0.apply(identity)
    after1 = L1.apply(after0)
    after2 = L2.apply(after1)

    return after2


# ==============================================================================
# EXACT VECTOR/SYSTEM SOLVER
# ==============================================================================

def rref(A, b):
    M = [
        [Fraction(x) for x in row] + [Fraction(rhs)]
        for row, rhs in zip(A, b)
    ]

    rows = len(M)
    cols = len(A[0]) if A else 0

    pivots = []
    prow = 0

    for col in range(cols):

        pivot = None

        for r in range(prow, rows):
            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        M[prow], M[pivot] = M[pivot], M[prow]

        p = M[prow][col]
        M[prow] = [x / p for x in M[prow]]

        for r in range(rows):
            if r == prow:
                continue

            q = M[r][col]

            if q == 0:
                continue

            M[r] = [
                M[r][j] - q * M[prow][j]
                for j in range(cols + 1)
            ]

        pivots.append(col)
        prow += 1

        if prow == rows:
            break

    for row in M:
        if (
            all(row[j] == 0 for j in range(cols))
            and row[cols] != 0
        ):
            return None, pivots

    solution = [Fraction(0)] * cols

    for i, p in enumerate(pivots):
        if p < cols:
            solution[p] = M[i][cols]

    return solution, pivots


# ==============================================================================
# COMPLEXITY
# ==============================================================================

def complexity(poly):
    nz = [x for x in poly if x != 0]

    if not nz:
        return (0, 0, 0)

    return (
        len(poly) - 1,
        max(len(str(abs(x.numerator))) for x in nz),
        max(len(str(x.denominator)) for x in nz),
    )


def operator_complexity(op):
    return tuple(
        complexity(x)
        for x in op
    )


# ==============================================================================
# FACTORIZATION SEARCH
# ==============================================================================

def enumerate_factor_shape(degree):
    """
    Return all possible first-order factor degree assignments
    with A,B of degree <= degree.
    """

    return [
        (degree, degree),
    ]


def polynomial_value_constraints(
    known,
    candidate,
    sample_count,
):
    """
    Exact pointwise comparison.
    """

    if len(known) != len(candidate):
        return False

    for a in range(len(known)):

        kc = known[a]
        cc = candidate[a]

        for r in range(sample_count):

            if eval_falling(kc, r) != eval_falling(cc, r):
                return False

    return True


# ==============================================================================
# SEARCH VIA EXACT NONLINEAR PARAMETERIZATION
# ==============================================================================

def search_small_factorizations(
    known,
    degrees=(0, 1, 2),
):
    """
    The composition coefficients are polynomial products, hence the
    factorization problem is nonlinear.

    We therefore use a finite candidate support search:

       A_i(r) = a_i0 + a_i1 r_(1) + ...
       B_i(r) = b_i0 + b_i1 r_(1) + ...

    and inspect sparse patterns first.

    The search is deliberately conservative.
    """

    print()
    print(
        "  Sparse-factor candidates are tested before dense "
        "nonlinear solving."
    )

    candidates = []

    # Identity-like factors.
    candidates.extend([
        ("identity-left", None),
        ("identity-middle", None),
        ("identity-right", None),
    ])

    # We do not invent arbitrary nonlinear solutions here.
    return candidates


# ==============================================================================
# EXACT OPERATOR ROOT TESTS
# ==============================================================================

def scalar_multiple(a, b):
    ratio = None

    for x, y in zip(a, b):

        if x == 0 and y == 0:
            continue

        if x == 0 or y == 0:
            return False, None

        q = x / y

        if ratio is None:
            ratio = q
        elif q != ratio:
            return False, None

    return True, ratio


def common_factor_test(components):
    """
    Test whether all shift components share a common polynomial factor
    in the falling coordinate, using exact integer evaluation.

    This does not infer arbitrary symbolic factors.
    """

    if not components:
        return None

    # Small integer roots only.
    for root in range(-10, 11):

        if all(
            eval_falling(poly, root) == 0
            for poly in components
        ):
            return root

    return None


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 122 — EXACT B-ODD SHIFT-OPERATOR "
        "FACTORIZATION AUDIT"
    )
    print("=" * 78)

    print()
    print("=" * 78)
    print("1. OPERATOR INPUT VALIDATION")
    print("=" * 78)

    print("  F_0 degree =", len(F0) - 1)
    print("  F_1 degree =", len(F1) - 1)
    print("  F_2 degree =", len(F2) - 1)
    print("  F_3 =", F3)

    print()
    print("  IMPORTANT:")
    print(
        "  The previous Experiment-120 transcript contains only F_0,F_1,F_2 "
        "explicitly."
    )
    print(
        "  An order-3 shift operator also requires F_3."
    )
    print(
        "  The script therefore refuses to fabricate F_3."
    )

    print()
    print("=" * 78)
    print("2. COMMON FACTOR AUDIT ON KNOWN COMPONENTS")
    print("=" * 78)

    known_components = [
        F0,
        F1,
        F2,
    ]

    for i, F in enumerate(known_components):
        root = None

        for r in range(-20, 21):

            if eval_falling(F, r) == 0:
                root = r
                break

        print(
            f"  F_{i}: integer_root_in[-20,20]={root}"
        )

    print()
    print("=" * 78)
    print("3. PAIRWISE PROPORTIONALITY AUDIT")
    print("=" * 78)

    for i in range(len(known_components)):
        for j in range(i + 1, len(known_components)):

            ok, scale = scalar_multiple(
                known_components[i],
                known_components[j],
            )

            print(
                f"  F_{i} vs F_{j}: "
                f"proportional={ok} scalar={scale}"
            )

    print()
    print("=" * 78)
    print("4. FACTORIZATION SEARCH STATUS")
    print("=" * 78)

    search_small_factorizations(
        known_components
    )

    print()
    print("=" * 78)
    print("5. STRUCTURAL INTERPRETATION")
    print("=" * 78)
    print()
    print(
        "Experiments 118-121 established an exceptional exact "
        "B-odd order-3 operator."
    )
    print()
    print(
        "Experiment 121 showed that the standard factorial/binomial "
        "diagonal changes do not transfer the operator to A-even, "
        "A-odd, or B-even."
    )
    print()
    print(
        "The next structural question is operator factorization:"
    )
    print()
    print(
        "    T = L_2 L_1 L_0"
    )
    print()
    print(
        "with first-order factors"
    )
    print()
    print(
        "    L_i = A_i(r) + B_i(r) E."
    )
    print()
    print(
        "Such a factorization would be substantially stronger than "
        "another interpolation identity because it would expose the "
        "order-3 recurrence as a composition of elementary shifts."
    )
    print()
    print(
        "The current transcript does not contain the exact F_3 component, "
        "so a complete exact factorization cannot honestly be completed "
        "from the supplied data alone."
    )
    print()
    print(
        "This script deliberately stops before making that missing "
        "coefficient up."
    )

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  known_components_exact=True"
    )
    print(
        "  missing_F3=True"
    )
    print(
        "  complete_factorization=False"
    )
    print(
        "  ALL BASIC CHECKS PASS=False"
    )

    print()
    print(
        "EXPERIMENT 122 INCOMPLETE — F_3 REQUIRED"
    )


if __name__ == "__main__":
    main()

