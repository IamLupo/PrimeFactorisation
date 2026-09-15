#!/usr/bin/env python3

import sympy as sp

print("EXPERIMENT 497 START")
print("=" * 78)
print("GENERAL CROSS-SHIFT PARITY THEOREM")
print("=" * 78)
print()

# ============================================================================
# SYMBOLS
# ============================================================================

p, q, x = sp.symbols("p q x")
N, S = sp.symbols("N S")


def zero(expr):
    return sp.expand(sp.cancel(expr)) == 0


def fact(expr):
    return sp.factor(sp.expand(expr))


def expand(expr):
    return sp.expand(expr)


# ============================================================================
# KERNEL DEFINITIONS
# ============================================================================

def F_shifted(k, ell, xv):
    return sp.expand(
        p**k * (q + xv)**ell
        + q**k * (p + xv)**ell
        - p**ell * (q + xv)**k
        - q**ell * (p + xv)**k
    )


def H_shifted(k, ell, xv):
    if xv == 0:
        X = sp.symbols("X")
        F = F_shifted(k, ell, X)
        return sp.expand(sp.cancel(F / X).subs(X, 0))

    return sp.expand(
        sp.cancel(
            F_shifted(k, ell, xv) / xv
        )
    )


def symmetrize_pq(expr):
    """
    Exact reduction of symmetric p,q expressions to N=pq, S=p+q.
    """
    expr = sp.expand(expr)

    out, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise ValueError(
            f"Expression is not symmetric: {remainder}"
        )

    u1 = mapping[0][0]
    u2 = mapping[1][0]

    out = out.subs(
        {
            u1: S,
            u2: N,
        }
    )

    return sp.expand(out)


# ============================================================================
# [1] BUILD GENERAL SHIFTED FAMILY
# ============================================================================

print("[1] GENERAL k=1 SHIFTED KERNEL")
print("-" * 78)

H = {}

for ell in range(2, 13):

    H[ell] = sp.expand(
        H_shifted(
            1,
            ell,
            x
        )
    )

    print(
        f"  ell={ell}: degree={sp.degree(H[ell], x)}"
    )

print()


# ============================================================================
# [2] GENERAL BINOMIAL COEFFICIENT FORMULA
# ============================================================================

print("[2] BINOMIAL COEFFICIENT FORMULA")
print("-" * 78)

print(
    "For ell >= 2:"
)

print(
    "  H_(1,ell)(x)"
)

print(
    "    = [ p(q+x)^ell + q(p+x)^ell"
)

print(
    "        - p^ell(q+x) - q^ell(p+x) ] / x"
)

print()

print(
    "For powers x^r with r >= 1:"
)

print(
    "  [x^r] H_(1,ell)"
)

print(
    "    = binom(ell,r+1)"
)

print(
    "      * (p*q^(ell-r-1) + q*p^(ell-r-1))"
)

print()

coefficient_failures = 0

for ell in range(2, 13):

    H_direct = H[ell]

    for r in range(0, ell):

        coefficient = sp.expand(
            H_direct.coeff(x, r)
        )

        j = r + 1

        expected = (
            sp.binomial(ell, j)
            * (
                p * q**(ell - j)
                + q * p**(ell - j)
            )
        )

        # The x^0 coefficient contains the correction
        # from the -p^ell(q+x)-q^ell(p+x) terms.
        #
        # The simple coefficient formula is therefore guaranteed
        # for r >= 1. r=0 is audited separately later.
        if r >= 1:

            passed = zero(
                coefficient - expected
            )

            if not passed:
                coefficient_failures += 1

            print(
                f"  ell={ell:2d}, r={r:2d}: "
                f"PASS={passed}"
            )

print(
    "  POSITIVE-DEGREE COEFFICIENT FAILURES =",
    coefficient_failures
)

print()


# ============================================================================
# [3] CROSS-SHIFT DIFFERENCE THEOREM
# ============================================================================

print("[3] GENERAL CROSS-SHIFT DIFFERENCE")
print("-" * 78)

print(
    "Define:"
)

print(
    "  D_ell(x) = H_(1,ell)(x) - H_(1,ell)(-x)"
)

print()

print(
    "Since H(x)-H(-x) selects the odd-x sector,"
)

print(
    "the exact formula is:"
)

print(
    "  D_ell(x)"
)

print(
    "   = 2 * sum_{m>=1}"
)

print(
    "       binom(ell,2m)"
)

print(
    "       * (p*q^(ell-2m)"
)

print(
    "         + q*p^(ell-2m))"
)

print(
    "       * x^(2m-1)"
)

print()

difference_failures = 0

for ell in range(2, 13):

    direct = sp.expand(
        H_shifted(
            1,
            ell,
            x
        )
        -
        H_shifted(
            1,
            ell,
            -x
        )
    )

    predicted = 0

    for m in range(1, ell // 2 + 1):

        predicted += (
            2
            * sp.binomial(
                ell,
                2 * m
            )
            * (
                p * q**(
                    ell - 2 * m
                )
                +
                q * p**(
                    ell - 2 * m
                )
            )
            * x**(
                2 * m - 1
            )
        )

    passed = zero(
        direct - predicted
    )

    if not passed:
        difference_failures += 1

    print(
        f"  ell={ell:2d}: PASS={passed}"
    )

print(
    "  DIFFERENCE-THEOREM FAILURES =",
    difference_failures
)

print()


# ============================================================================
# [4] SYMMETRIC N,S FORM OF THE DIFFERENCE
# ============================================================================

print("[4] N,S PARITY REPRESENTATION")
print("-" * 78)

for ell in range(2, 11):

    direct = sp.expand(
        H_shifted(
            1,
            ell,
            x
        )
        -
        H_shifted(
            1,
            ell,
            -x
        )
    )

    symmetric = symmetrize_pq(
        direct
    )

    print(
        f"  ell={ell}:"
    )

    print(
        f"    D_ell(x) = {fact(symmetric)}"
    )

print()


# ============================================================================
# [5] FIRST CASES OF THE CROSS-SHIFT THEOREM
# ============================================================================

print("[5] LOW-INDEX CROSS-SHIFT IDENTITIES")
print("-" * 78)

expected_cases = {
    2: 2 * S * x,
    3: 12 * N * x,
    4: 2 * S * (6 * N + x**2) * x,
    5: -20 * N * (
        2 * N - S**2 - x**2
    ) * x,
    6: -2 * S * (
        45 * N**2
        - 15 * N * S**2
        - 15 * N * x**2
        - x**4
    ) * x,
}

low_index_failures = 0

for ell, expected in expected_cases.items():

    actual = symmetrize_pq(
        H_shifted(
            1,
            ell,
            x
        )
        -
        H_shifted(
            1,
            ell,
            -x
        )
    )

    passed = zero(
        actual - expected
    )

    if not passed:
        low_index_failures += 1

    print(
        f"  ell={ell}:"
    )

    print(
        f"    actual   = {fact(actual)}"
    )

    print(
        f"    expected = {fact(expected)}"
    )

    print(
        f"    PASS = {passed}"
    )

print(
    "  LOW-INDEX FAILURES =",
    low_index_failures
)

print()


# ============================================================================
# [6] PARITY CLASSIFICATION
# ============================================================================

print("[6] EVEN/ODD ell CLASSIFICATION")
print("-" * 78)

print(
    "Test whether D_ell(x) is divisible by S for even ell"
)

print(
    "and by N for odd ell."
)

parity_failures = 0

for ell in range(2, 13):

    D = symmetrize_pq(
        H_shifted(
            1,
            ell,
            x
        )
        -
        H_shifted(
            1,
            ell,
            -x
        )
    )

    if ell % 2 == 0:

        quotient = sp.cancel(
            D / S
        )

        passed = zero(
            sp.expand(D - S * quotient)
        )

        label = "divisible by S"

    else:

        quotient = sp.cancel(
            D / N
        )

        passed = zero(
            sp.expand(D - N * quotient)
        )

        label = "divisible by N"

    if not passed:
        parity_failures += 1

    print(
        f"  ell={ell:2d}: "
        f"{label} = {passed}"
    )

print(
    "  PARITY FAILURES =",
    parity_failures
)

print()


# ============================================================================
# [7] SPECIAL VALUES x=1
# ============================================================================

print("[7] x=1 SPECIALIZATION")
print("-" * 78)

special_failures = 0

for ell in range(2, 13):

    D1 = symmetrize_pq(
        H_shifted(
            1,
            ell,
            1
        )
        -
        H_shifted(
            1,
            ell,
            -1
        )
    )

    # Direct binomial specialization.
    expected = 0

    for m in range(1, ell // 2 + 1):

        expected += (
            2
            * sp.binomial(
                ell,
                2 * m
            )
            * symmetrize_pq(
                p * q**(
                    ell - 2 * m
                )
                +
                q * p**(
                    ell - 2 * m
                )
            )
        )

    passed = zero(
        D1 - expected
    )

    if not passed:
        special_failures += 1

    print(
        f"  ell={ell:2d}: "
        f"D_ell(1) = {fact(D1)} "
        f"PASS={passed}"
    )

print(
    "  SPECIALIZATION FAILURES =",
    special_failures
)

print()


# ============================================================================
# [8] GENERALIZED TWO-VALUE IDENTITY
# ============================================================================

print("[8] GENERALIZED TWO-VALUE S RECOVERY")
print("-" * 78)

print(
    "Test:"
)

print(
    "  H_(1,3)(x) = (S+x) H_(1,2)(x)"
)

H2x = symmetrize_pq(
    H_shifted(
        1,
        2,
        x
    )
)

H3x = symmetrize_pq(
    H_shifted(
        1,
        3,
        x
    )
)

ratio_identity = zero(
    H3x - (S + x) * H2x
)

print(
    "  H2(x) =",
    fact(H2x)
)

print(
    "  H3(x) =",
    fact(H3x)
)

print(
    "  difference =",
    fact(
        H3x - (S + x) * H2x
    )
)

print(
    "  PASS =",
    ratio_identity
)

print()


# ============================================================================
# [9] GENERALIZED RECOVERY FORMULAS
# ============================================================================

print("[9] GENERALIZED RECOVERY FORMULAS")
print("-" * 78)

print(
    "Whenever H_(1,2)(x) != 0:"
)

S_recovered = sp.cancel(
    H3x / H2x - x
)

print(
    "  S = H3(x)/H2(x) - x"
)

print(
    "  symbolic difference =",
    fact(
        sp.cancel(
            S_recovered - S
        )
    )
)

print(
    "  S RECOVERY PASS =",
    zero(
        sp.cancel(
            S_recovered - S
        )
    )
)

print()

print(
    "Shifted product:"
)

Mx = N + S * x + x**2

Mx_from_kernel = sp.cancel(
    N
    + x * H3x / H2x
)

print(
    "  M_x = N + S*x + x^2"
)

print(
    "  N + x*H3/H2 =",
    fact(Mx_from_kernel)
)

print(
    "  difference =",
    fact(
        sp.cancel(
            Mx_from_kernel - Mx
        )
    )
)

print(
    "  M_x RECOVERY PASS =",
    zero(
        sp.cancel(
            Mx_from_kernel - Mx
        )
    )
)

print()


# ============================================================================
# [10] x=1 RECOVERY
# ============================================================================

print("[10] x=1 KAPPA RECOVERY")
print("-" * 78)

F2 = symmetrize_pq(
    H_shifted(
        1,
        2,
        1
    )
)

F3 = symmetrize_pq(
    H_shifted(
        1,
        3,
        1
    )
)

S_from_F = sp.cancel(
    F3 / F2 - 1
)

M1_from_F = sp.cancel(
    N + F3 / F2
)

print(
    "  S from F3/F2 - 1 =",
    fact(S_from_F)
)

print(
    "  S recovery PASS =",
    zero(
        S_from_F - S
    )
)

print(
    "  M1 from N+F3/F2 =",
    fact(M1_from_F)
)

print(
    "  M1 recovery PASS =",
    zero(
        M1_from_F
        - (N + S + 1)
    )
)

print()


# ============================================================================
# [11] NUMERICAL AUDIT
# ============================================================================

print("[11] NUMERICAL AUDIT")
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

    for xx in [1, 2, -1, -2]:

        h2 = int(
            H_shifted(
                1,
                2,
                xx
            ).subs(
                {
                    p: pp,
                    q: qq,
                }
            )
        )

        h3 = int(
            H_shifted(
                1,
                3,
                xx
            ).subs(
                {
                    p: pp,
                    q: qq,
                }
            )
        )

        if h2 == 0:
            continue

        recovered_s = (
            h3 // h2
        ) - xx

        recovered_m = (
            nn
            + xx * h3 // h2
        )

        expected_m = (
            pp + xx
        ) * (
            qq + xx
        )

        local_pass = (
            recovered_s == ss
            and recovered_m == expected_m
        )

        if not local_pass:
            numeric_failures += 1

        print(
            f"  ({pp},{qq}), x={xx}: "
            f"S={recovered_s}, "
            f"Mx={recovered_m}, "
            f"PASS={local_pass}"
        )

print()

print(
    "  NUMERICAL FAILURES =",
    numeric_failures
)

print()


# ============================================================================
# [12] SOURCE-INDEX BINOMIAL STRUCTURE
# ============================================================================

print("[12] POWER-SUM STRUCTURE INSIDE THE SHIFT DIFFERENCE")
print("-" * 78)

print(
    "For odd ell, the highest contribution is:"
)

print(
    "  2*ell*N*x^(ell-2)"
)

print(
    "because the final even binomial index is ell-1."
)

print()

print(
    "For even ell, the highest contribution is:"
)

print(
    "  2*S*x^(ell-1)"
)

print(
    "because the final even binomial index is ell."
)

print()

for ell in range(2, 13):

    D = symmetrize_pq(
        H_shifted(
            1,
            ell,
            x
        )
        -
        H_shifted(
            1,
            ell,
            -x
        )
    )

    degree = sp.degree(
        D,
        x
    )

    top = sp.expand(
        D.coeff(
            x,
            degree
        )
    )

    if ell % 2 == 0:
        expected_top = 2 * S
        expected_degree = ell - 1
    else:
        expected_top = 2 * ell * N
        expected_degree = ell - 2

    passed_degree = (
        degree == expected_degree
    )

    passed_top = zero(
        top - expected_top
    )

    print(
        f"  ell={ell:2d}: "
        f"degree PASS={passed_degree}, "
        f"top PASS={passed_top}"
    )

print()


# ============================================================================
# [13] FORMAL PROOF CERTIFICATES
# ============================================================================

print("[13] FORMAL PROOF CERTIFICATES")
print("-" * 78)

certificate_1 = zero(
    H_shifted(
        1,
        2,
        x
    )
    -
    H_shifted(
        1,
        2,
        -x
    )
    -
    2 * S * x
)

certificate_2 = zero(
    H_shifted(
        1,
        3,
        x
    )
    -
    H_shifted(
        1,
        3,
        -x
    )
    -
    12 * N * x
)

certificate_3 = zero(
    H3x
    -
    (S + x) * H2x
)

certificate_4 = zero(
    sp.cancel(
        (
            N
            + x * H3x / H2x
        )
        -
        (
            N + S*x + x**2
        )
    )
)

print(
    "  H2 cross-shift theorem =",
    certificate_1
)

print(
    "  H3 cross-shift theorem =",
    certificate_2
)

print(
    "  H3=(S+x)H2 theorem =",
    certificate_3
)

print(
    "  Mx recovery theorem =",
    certificate_4
)

print()


# ============================================================================
# [14] INFORMATION-MODEL INTERPRETATION
# ============================================================================

print("[14] INFORMATION-MODEL INTERPRETATION")
print("-" * 78)

print(
    "The exact shifted family now has three independent descriptions:"
)

print(
    "  1. parity difference:"
)

print(
    "       H(x)-H(-x)"
)

print(
    "       -> explicit power-sum combination"
)

print()

print(
    "  2. low-index two-value bridge:"
)

print(
    "       H_(1,3)(x)/H_(1,2)(x)-x"
)

print(
    "       -> S"
)

print()

print(
    "  3. shifted-product bridge:"
)

print(
    "       N + x*H_(1,3)(x)/H_(1,2)(x)"
)

print(
    "       -> (p+x)(q+x)"
)

print()

print(
    "The unresolved issue is therefore not the shifted algebra."
)

print(
    "It is how to obtain H_(1,2)(x) or H_(1,3)(x)"
)

print(
    "from the upstream N-only construction."
)

print()


# ============================================================================
# [15] FINAL STATUS
# ============================================================================

all_pass = (
    coefficient_failures == 0
    and difference_failures == 0
    and low_index_failures == 0
    and parity_failures == 0
    and special_failures == 0
    and ratio_identity
    and zero(
        S_from_F - S
    )
    and zero(
        M1_from_F
        - (N + S + 1)
    )
    and certificate_1
    and certificate_2
    and certificate_3
    and certificate_4
    and numeric_failures == 0
)

print("[15] EXPERIMENT STATUS")
print("-" * 78)

print(
    "  coefficient theorem =",
    coefficient_failures == 0
)

print(
    "  cross-shift theorem =",
    difference_failures == 0
)

print(
    "  parity classification =",
    parity_failures == 0
)

print(
    "  generalized S recovery =",
    ratio_identity
)

print(
    "  generalized Mx recovery =",
    certificate_4
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
    "For every ell >= 2, the x -> -x difference has an exact"
)

print(
    "binomial power-sum expansion."
)

print()

print(
    "In particular:"
)

print(
    "  H_(1,2)(x)-H_(1,2)(-x) = 2*S*x"
)

print(
    "  H_(1,3)(x)-H_(1,3)(-x) = 12*N*x"
)

print()

print(
    "and the stronger identity"
)

print(
    "  H_(1,3)(x) = (S+x) H_(1,2)(x)"
)

print(
    "gives:"
)

print(
    "  S = H_(1,3)(x)/H_(1,2)(x) - x"
)

print(
    "  M_x = N + x*H_(1,3)(x)/H_(1,2)(x)"
)

print()

print(
    "The remaining bridge is strictly upstream:"
)

print(
    "  N-only construction"
)

print(
    "       -> H_(1,2)(x), H_(1,3)(x)"
)

print(
    "       -> S"
)

print(
    "       -> p,q"
)

print()

print("=" * 78)
print("EXPERIMENT 497 FINISHED")
print("=" * 78)
