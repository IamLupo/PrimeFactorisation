#!/usr/bin/env python3

import sympy as sp

print("EXPERIMENT 498 START")
print("=" * 78)
print("SHIFTED-KERNEL FIELD CLOSURE AND S-ELIMINATION AUDIT")
print("=" * 78)
print()

# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")
x, a, b = sp.symbols("x a b")

Y = sp.symbols("Y")


def zero(expr):
    return sp.expand(sp.cancel(expr)) == 0


def fact(expr):
    return sp.factor(sp.expand(expr))


def symmetrize_pq(expr):
    expr = sp.expand(expr)

    sym_expr, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise ValueError(
            f"Expression is not symmetric in p,q: {remainder}"
        )

    u1 = mapping[0][0]
    u2 = mapping[1][0]

    return sp.expand(
        sym_expr.subs(
            {
                u1: S,
                u2: N,
            }
        )
    )


# ============================================================================
# SHIFTED KERNEL
# ============================================================================

def F_shifted(k, ell, xv):
    return sp.expand(
        p**k * (q + xv)**ell
        + q**k * (p + xv)**ell
        - p**ell * (q + xv)**k
        - q**ell * (p + xv)**k
    )


def H12(xv):
    return symmetrize_pq(
        sp.cancel(
            F_shifted(1, 2, xv) / xv
        )
    )


def H13(xv):
    return symmetrize_pq(
        sp.cancel(
            F_shifted(1, 3, xv) / xv
        )
    )


# ============================================================================
# [1] BASIC GENERATORS
# ============================================================================

print("[1] BASIC SHIFTED GENERATORS")
print("-" * 78)

h2x = sp.expand(H12(x))
h3x = sp.expand(H13(x))

print("  H2(x) =", fact(h2x))
print("  H3(x) =", fact(h3x))

print()

print(
    "  H2(x) = 6*N - S^2 + S*x"
)

print(
    "  H3(x) = (S+x)*H2(x)"
)

print()


# ============================================================================
# [2] GENERAL TWO-SHIFT DIFFERENCE
# ============================================================================

print("[2] TWO-SHIFT DIFFERENCE THEOREM")
print("-" * 78)

h2a = sp.expand(H12(a))
h2b = sp.expand(H12(b))

difference = sp.expand(h2a - h2b)
expected_difference = S * (a - b)

print("  H2(a) =", fact(h2a))
print("  H2(b) =", fact(h2b))
print("  H2(a)-H2(b) =", fact(difference))
print("  Expected      =", fact(expected_difference))
print(
    "  PASS =",
    zero(difference - expected_difference)
)

print()


# ============================================================================
# [3] GENERAL TWO-SHIFT RECOVERY
# ============================================================================

print("[3] GENERAL TWO-SHIFT RECOVERY")
print("-" * 78)

S_two_shift = sp.cancel(
    (h2a - h2b) / (a - b)
)

print(
    "  S = (H2(a)-H2(b))/(a-b)"
)

print(
    "  recovered S =",
    fact(S_two_shift)
)

print(
    "  difference =",
    fact(
        sp.cancel(
            S_two_shift - S
        )
    )
)

print(
    "  PASS =",
    zero(
        sp.cancel(
            S_two_shift - S
        )
    )
)

print()


# ============================================================================
# [4] PRODUCT OF TWO SHIFTED H2 VALUES
# ============================================================================

print("[4] PRODUCT CANCELLATION SEARCH")
print("-" * 78)

product_ab = sp.expand(
    h2a * h2b
)

print(
    "  H2(a)*H2(b) =",
    fact(product_ab)
)

print()

print(
    "Search for combinations using:"
)

print(
    "  H2(a)+H2(b)"
)

print(
    "  H2(a)-H2(b)"
)

print(
    "  H2(a)*H2(b)"
)

print(
    "  H2(a)^2 + H2(b)^2"
)

print(
    "  H2(a)^2 - H2(b)^2"
)

print()

combos = {
    "sum":
        sp.expand(h2a + h2b),

    "difference":
        sp.expand(h2a - h2b),

    "product":
        sp.expand(h2a * h2b),

    "square_sum":
        sp.expand(h2a**2 + h2b**2),

    "square_difference":
        sp.expand(h2a**2 - h2b**2),
}

for name, expr in combos.items():

    contains_S = (
        sp.diff(expr, S) != 0
    )

    print(
        f"  {name:18s}: "
        f"S-dependent={contains_S}"
    )

print()


# ============================================================================
# [5] DISCRIMINANT OF THE SHIFTED H2 QUADRATIC
# ============================================================================

print("[5] SHIFTED-EVALUATION QUADRATIC")
print("-" * 78)

print(
    "Treat H2(x) as a linear polynomial in x:"
)

print(
    "  H2(x) = S*x + (6N-S^2)"
)

poly_h2 = sp.Poly(
    h2x,
    x,
)

print(
    "  coefficients =",
    poly_h2.all_coeffs()
)

print()

intercept = poly_h2.coeff_monomial(
    x**0
)

slope = poly_h2.coeff_monomial(
    x
)

print(
    "  slope     =", slope
)

print(
    "  intercept =", intercept
)

print(
    "  slope - S =",
    fact(slope - S)
)

print(
    "  intercept + S^2 =",
    fact(intercept + S**2)
)

print()


# ============================================================================
# [6] CAN TWO H2 VALUES GIVE A PURE-N QUANTITY WITHOUT S?
# ============================================================================

print("[6] ALGEBRAIC S-ELIMINATION FROM TWO H2 VALUES")
print("-" * 78)

u, v = sp.symbols("u v")

eq_a = sp.expand(
    u - h2a
)

eq_b = sp.expand(
    v - h2b
)

print(
    "  equations:"
)

print(
    "    u =", fact(h2a)
)

print(
    "    v =", fact(h2b)
)

print()

resultant_uv = sp.factor(
    sp.resultant(
        eq_a,
        eq_b,
        S
    )
)

print(
    "  resultant eliminating S ="
)

print(
    "   ",
    resultant_uv
)

print()

print(
    "  The resultant is a relation among N,u,v."
)

print()


# ============================================================================
# [7] SPECIAL SHIFTS a=1,b=-1
# ============================================================================

print("[7] SPECIAL SHIFTS +1 / -1")
print("-" * 78)

hp = sp.expand(
    H12(1)
)

hm = sp.expand(
    H12(-1)
)

print(
    "  H(+1) =", fact(hp)
)

print(
    "  H(-1) =", fact(hm)
)

print(
    "  sum   =", fact(hp + hm)
)

print(
    "  diff  =", fact(hp - hm)
)

print(
    "  diff/2 =", fact((hp - hm) / 2)
)

print()


# ============================================================================
# [8] DISCRIMINANT-LIKE PRODUCT FOR +1/-1
# ============================================================================

print("[8] PRODUCT OF H(+1),H(-1)")
print("-" * 78)

product_pm = sp.expand(
    hp * hm
)

print(
    "  H(+1)H(-1) =",
    fact(product_pm)
)

print()

print(
    "  Compare against candidate:"
)

candidates = {
    "N^2":
        N**2,

    "N*(N+1)":
        N * (N + 1),

    "N^2 - 6N":
        N**2 - 6*N,

    "(6N-S^2)^2-S^2":
        (6*N - S**2)**2 - S**2,

    "(6N-S^2)^2":
        (6*N - S**2)**2,
}

for name, candidate in candidates.items():

    difference_candidate = sp.expand(
        product_pm - candidate
    )

    print(
        f"  {name:25s}: "
        f"difference = {fact(difference_candidate)}"
    )

print()


# ============================================================================
# [9] F2 / H(-1) ALGEBRA
# ============================================================================

print("[9] x=1 KAPPA VS x=-1 TRANSFER ALGEBRA")
print("-" * 78)

F2 = hp
Hminus2 = hm

F3 = symmetrize_pq(
    sp.cancel(
        F_shifted(
            1,
            3,
            1
        )
    )
)

Hminus3 = symmetrize_pq(
    sp.cancel(
        F_shifted(
            1,
            3,
            -1
        )
        / (-1)
    )
)

print(
    "  F2       =", fact(F2)
)

print(
    "  F3       =", fact(F3)
)

print(
    "  Hminus2  =", fact(Hminus2)
)

print(
    "  Hminus3  =", fact(Hminus3)
)

print()

relation_2 = sp.expand(
    F2 - Hminus2
)

relation_3 = sp.expand(
    F3 - Hminus3
)

print(
    "  F2-Hminus2 =",
    fact(relation_2)
)

print(
    "  F3-Hminus3 =",
    fact(relation_3)
)

print()


# ============================================================================
# [10] SEARCH FOR N-ONLY TRANSFER POLYNOMIALS
# ============================================================================

print("[10] N-ONLY TRANSFER POLYNOMIAL SEARCH")
print("-" * 78)

print(
    "Search whether Hminus2 or Hminus3 can be represented as"
)

print(
    "a polynomial in N multiplied by F2/F3 combinations."
)

basis = []

for max_degree in [0, 1, 2, 3]:

    current = []

    for i in range(2, 4):

        Fi = {
            2: F2,
            3: F3,
        }[i]

        for d in range(max_degree + 1):

            current.append(
                (f"N^{d}*F{i}", N**d * Fi)
            )

    basis = current

    print(
        f"  coefficient degree <= {max_degree}:"
    )

    # Solve exact polynomial linear combination
    # against target Hminus2.
    coeffs = sp.symbols(
        f"c0:{len(basis)}"
    )

    candidate = sum(
        c * expr
        for c, (_, expr)
        in zip(
            coeffs,
            basis
        )
    )

    polynomial = sp.Poly(
        sp.expand(candidate - Hminus2),
        N,
        S
    )

    equations = [
        coefficient
        for monomial, coefficient
        in polynomial.terms()
    ]

    if equations:

        sol = sp.linsolve(
            equations,
            coeffs
        )

        hits = []

        for solution in sol:

            if any(
                value != 0
                for value in solution
            ):

                hits.append(
                    solution
                )

        print(
            "    Hminus2 transfer hits =",
            len(hits)
        )

    else:

        print(
            "    Hminus2 transfer hits = 0"
        )

print()


# ============================================================================
# [11] GENERIC SHIFTED-PRODUCT IDENTITY
# ============================================================================

print("[11] SHIFTED PRODUCT RECOVERY")
print("-" * 78)

Mx = N + S*x + x**2

Mx_from_two_shift = sp.expand(
    N
    + x * (
        H13(x)
        / H12(x)
    )
)

print(
    "  M_x =", fact(Mx)
)

print(
    "  reconstructed =",
    fact(
        sp.cancel(
            Mx_from_two_shift
        )
    )
)

print(
    "  difference =",
    fact(
        sp.cancel(
            Mx_from_two_shift - Mx
        )
    )
)

print(
    "  PASS =",
    zero(
        sp.cancel(
            Mx_from_two_shift - Mx
        )
    )
)

print()


# ============================================================================
# [12] H3/H2 STRUCTURE FOR ARBITRARY x
# ============================================================================

print("[12] RATIO STRUCTURE")
print("-" * 78)

ratio_x = sp.cancel(
    H13(x) / H12(x)
)

print(
    "  H3(x)/H2(x) =",
    fact(ratio_x)
)

print(
    "  ratio - (S+x) =",
    fact(
        sp.cancel(
            ratio_x - (S+x)
        )
    )
)

print()


# ============================================================================
# [13] NUMERICAL AUDIT
# ============================================================================

print("[13] NUMERICAL AUDIT")
print("-" * 78)

instances = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),
]

numeric_failures = 0

for pp, qq in instances:

    nn = pp * qq
    ss = pp + qq

    for xx in [1, -1, 2, -2, 3]:

        h2v = int(
            H12(xx).subs(
                {
                    N: nn,
                    S: ss,
                }
            )
        )

        h3v = int(
            H13(xx).subs(
                {
                    N: nn,
                    S: ss,
                }
            )
        )

        if h2v == 0:
            print(
                f"  ({pp},{qq}), x={xx}: "
                f"H2=0 SKIPPED"
            )
            continue

        recovered_s = (
            h3v // h2v
        ) - xx

        recovered_m = (
            nn
            + xx * h3v // h2v
        )

        expected_m = (
            pp + xx
        ) * (
            qq + xx
        )

        passed = (
            recovered_s == ss
            and recovered_m == expected_m
        )

        if not passed:
            numeric_failures += 1

        print(
            f"  ({pp},{qq}), x={xx}: "
            f"S={recovered_s}, "
            f"Mx={recovered_m}, "
            f"PASS={passed}"
        )

print()

print(
    "  NUMERICAL FAILURES =",
    numeric_failures
)

print()


# ============================================================================
# [14] STRUCTURAL CONCLUSION
# ============================================================================

print("[14] STRUCTURAL CONCLUSION")
print("-" * 78)

print(
    "The shifted H2 family is exactly:"
)

print(
    "  H2(x)=6N-S^2+Sx."
)

print()

print(
    "Therefore its entire information content is"
)

print(
    "equivalent to the pair"
)

print(
    "  (S, 6N-S^2)."
)

print()

print(
    "Two shifted evaluations recover S directly:"
)

print(
    "  S=(H2(a)-H2(b))/(a-b)."
)

print()

print(
    "The H3 value adds no independent hidden variable because:"
)

print(
    "  H3(x)=(S+x)H2(x)."
)

print()

print(
    "Hence the unresolved upstream problem can be stated even more sharply:"
)

print(
    "  N-only construction"
)

print(
    "       -> ONE nontrivial source-free scalar"
)

print(
    "       -> determines H2(x)"
)

print(
    "       -> S"
)

print(
    "       -> p,q"
)

print()

print(
    "Equivalently, it is enough to generate"
)

print(
    "  T(N,p,q)=6N-S^2+Sx"
)

print(
    "for one fixed nonzero x."
)

print()


# ============================================================================
# [15] STATUS
# ============================================================================

ratio_pass = zero(
    sp.cancel(
        ratio_x - (S + x)
    )
)

mx_pass = zero(
    sp.cancel(
        Mx_from_two_shift - Mx
    )
)

two_shift_pass = zero(
    S_two_shift - S
)

all_pass = (
    zero(
        difference
        - expected_difference
    )
    and two_shift_pass
    and ratio_pass
    and mx_pass
    and numeric_failures == 0
)

print("[15] EXPERIMENT STATUS")
print("-" * 78)

print(
    "  two-shift difference theorem =",
    zero(
        difference
        - expected_difference
    )
)

print(
    "  two-shift S recovery =",
    two_shift_pass
)

print(
    "  H3/H2 ratio theorem =",
    ratio_pass
)

print(
    "  shifted-product recovery =",
    mx_pass
)

print(
    "  numerical audit =",
    numeric_failures == 0
)

print(
    "  OVERALL EXACT AUDIT =",
    all_pass
)

print()

print(
    "MAIN RESULT"
)

print(
    "-" * 78
)

print(
    "Experiment 498 reduces the shifted construction to a minimal"
)

print(
    "information interface."
)

print()

print(
    "For H2(x)=6N-S^2+Sx:"
)

print(
    "  H2(a)-H2(b)=S(a-b)."
)

print(
    "Therefore two shifted evaluations recover S exactly."
)

print()

print(
    "Furthermore:"
)

print(
    "  H3(x)=(S+x)H2(x),"
)

print(
    "so H3 introduces no new independent information."
)

print()

print(
    "The next upstream target is consequently:"
)

print(
    "  derive H2(x) or an algebraically equivalent scalar"
)

print(
    "  from the existing N-only homogeneous/KAPPA construction."
)

print()

print(
    "This experiment deliberately does NOT perform large"
)

print(
    "resultant searches or high-degree determinant expansion."
)

print(
    "The remaining problem is now an interface-identification"
)

print(
    "problem rather than another downstream algebra problem."
)

print()
print("=" * 78)
print("EXPERIMENT 498 FINISHED")
print("=" * 78)
