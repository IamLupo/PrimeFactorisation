#!/usr/bin/env python3

import sympy as sp


print("EXPERIMENT 492 START")
print("=" * 78)
print("SHIFTED-EVALUATION SYMMETRY AND N-ONLY INVARIANT SEARCH")
print("=" * 78)
print()


# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
N, S, x = sp.symbols("N S x")


# ============================================================================
# HELPERS
# ============================================================================

def zero(expr):
    return sp.expand(expr) == 0


def fact(expr):
    return sp.factor(sp.expand(expr))


def cancel_fact(expr):
    return sp.factor(
        sp.cancel(
            sp.together(expr)
        )
    )


# ============================================================================
# SYMMETRIC REDUCTION
# ============================================================================

def symmetric_reduce(expr):
    """
    Reduce a symmetric polynomial in p,q to N=pq and S=p+q.

    We use q=S-p and reduce modulo p^2-S*p+N.
    """
    expr_sub = sp.expand(
        expr.subs(q, S - p)
    )

    poly = sp.Poly(
        expr_sub,
        p,
        domain="EX"
    )

    modulus = sp.Poly(
        p**2 - S*p + N,
        p,
        domain="EX"
    )

    rem = sp.rem(
        poly,
        modulus
    ).as_expr()

    return sp.factor(
        sp.expand(rem)
    )


# ============================================================================
# [1] GENERAL SHIFTED KERNEL
# ============================================================================

print("[1] GENERAL SHIFTED KERNEL")
print("-" * 78)

def F_shifted(k, ell, xv):
    return sp.expand(
        p**k * (q + xv)**ell
        + q**k * (p + xv)**ell
        - p**ell * (q + xv)**k
        - q**ell * (p + xv)**k
    )


def H_shifted(k, ell, xv):
    """
    H = F/x.
    The x=0 value is defined by exact polynomial division.
    """

    F = sp.Poly(
        F_shifted(k, ell, x),
        x,
        domain="EX"
    )

    quotient, remainder = sp.div(
        F,
        sp.Poly(x, x)
    )

    if not remainder.is_zero:
        raise ValueError(
            "Expected F to be divisible by x."
        )

    H = quotient.as_expr()

    if xv == x:
        return sp.expand(H)

    return sp.expand(
        H.subs(x, xv)
    )


# ============================================================================
# [2] BASIC H_(1,l) POLYNOMIALS
# ============================================================================

print("[2] k=1 SHIFTED POLYNOMIALS")
print("-" * 78)

ell_values = [2, 3, 4, 5, 6]

H_polys = {}

for ell in ell_values:

    Hpq = H_shifted(
        1,
        ell,
        x
    )

    Hns = symmetric_reduce(
        Hpq
    )

    H_polys[ell] = Hns

    print(
        f"  ell={ell}:"
    )

    print(
        f"    H(x) = {Hns}"
    )

print()


# ============================================================================
# [3] PARITY DECOMPOSITION IN x
# ============================================================================

print("[3] x-PARITY DECOMPOSITION")
print("-" * 78)

parity_data = {}

for ell, Hexpr in H_polys.items():

    even_part = sp.expand(
        (Hexpr + Hexpr.subs(x, -x)) / 2
    )

    odd_part = sp.expand(
        (Hexpr - Hexpr.subs(x, -x)) / (2*x)
    )

    even_part = fact(even_part)
    odd_part = fact(odd_part)

    parity_data[ell] = (
        even_part,
        odd_part
    )

    print(
        f"  ell={ell}:"
    )

    print(
        f"    even H = {even_part}"
    )

    print(
        f"    odd/x  = {odd_part}"
    )

print()


# ============================================================================
# [4] EXACT x = +/- a EVALUATIONS
# ============================================================================

print("[4] PAIRED SHIFT EVALUATIONS")
print("-" * 78)

shifts = [-2, -1, 1, 2]

paired = {}

for ell in ell_values:

    print(
        f"  ell={ell}:"
    )

    for a in [1, 2]:

        H_plus = fact(
            H_polys[ell].subs(
                x,
                a
            )
        )

        H_minus = fact(
            H_polys[ell].subs(
                x,
                -a
            )
        )

        pair_sum = fact(
            H_plus + H_minus
        )

        pair_diff = fact(
            H_plus - H_minus
        )

        paired[(ell, a)] = (
            H_plus,
            H_minus,
            pair_sum,
            pair_diff
        )

        print(
            f"    a={a}:"
        )

        print(
            f"      H(+a) = {H_plus}"
        )

        print(
            f"      H(-a) = {H_minus}"
        )

        print(
            f"      sum    = {pair_sum}"
        )

        print(
            f"      diff   = {pair_diff}"
        )

    print()

# ============================================================================
# [5] N-ONLY TEST OF PAIRED COMBINATIONS
# ============================================================================

print("[5] N-ONLY TEST OF PAIRED COMBINATIONS")
print("-" * 78)

pair_candidates = {}

for key, (
    Hplus,
    Hminus,
    pair_sum,
    pair_diff
) in paired.items():

    ell, a = key

    candidates = {
        f"H{ell}(+{a})+H{ell}(-{a})":
            pair_sum,

        f"H{ell}(+{a})-H{ell}(-{a})":
            pair_diff,

        f"(H+ + H-)/2":
            sp.expand(pair_sum / 2),

        f"(H+ - H-)/(2a)":
            sp.expand(pair_diff / (2*a)),
    }

    for name, expr in candidates.items():

        expr = fact(expr)

        derivative = sp.factor(
            sp.diff(
                expr,
                S
            )
        )

        pure_N = zero(
            derivative
        )

        pair_candidates[
            name
        ] = expr

        print(
            f"  {name}"
        )

        print(
            f"    expr    = {expr}"
        )

        print(
            f"    pure-N  = {pure_N}"
        )

    print()


# ============================================================================
# [6] SEARCH SIMPLE QUADRATIC COMBINATIONS
# ============================================================================

print("[6] QUADRATIC PAIRED-INVARIANT SEARCH")
print("-" * 78)

# Select the compactest pair features.
feature_names = []
feature_values = []

for ell in [2, 3, 4]:

    for a in [1, 2]:

        Hplus = paired[(ell, a)][0]
        Hminus = paired[(ell, a)][1]

        feature_names.extend(
            [
                f"H{ell}(+{a})",
                f"H{ell}(-{a})",
            ]
        )

        feature_values.extend(
            [
                Hplus,
                Hminus,
            ]
        )

quadratic_candidates = []

# Include simple products of feature pairs.
for i in range(len(feature_values)):

    for j in range(i, len(feature_values)):

        expr = sp.expand(
            feature_values[i]
            * feature_values[j]
        )

        derivative = sp.diff(
            expr,
            S
        )

        if zero(derivative):

            quadratic_candidates.append(
                (
                    f"{feature_names[i]}*{feature_names[j]}",
                    fact(expr)
                )
            )

print(
    "  pure-N quadratic products found =",
    len(quadratic_candidates)
)

for name, expr in quadratic_candidates[:30]:

    print(
        f"  {name}:"
    )

    print(
        f"    {expr}"
    )

print()


# ============================================================================
# [7] DISCRIMINANT-LIKE PAIRED COMBINATIONS
# ============================================================================

print("[7] DISCRIMINANT-LIKE COMBINATIONS")
print("-" * 78)

discriminant_candidates = []

for ell in [2, 3, 4, 5]:

    for a in [1, 2]:

        Hplus, Hminus, _, _ = paired[
            (ell, a)
        ]

        # Basic two-point symmetric invariants.
        candidates = {
            "sum^2 - diff^2":
                sp.expand(
                    (Hplus + Hminus)**2
                    - (Hplus - Hminus)**2
                ),

            "product":
                sp.expand(
                    Hplus * Hminus
                ),

            "sum":
                sp.expand(
                    Hplus + Hminus
                ),

            "difference":
                sp.expand(
                    Hplus - Hminus
                ),
        }

        for label, expr in candidates.items():

            expr = fact(expr)

            if zero(
                sp.diff(
                    expr,
                    S
                )
            ):

                discriminant_candidates.append(
                    (
                        ell,
                        a,
                        label,
                        expr
                    )
                )

print(
    "  pure-N discriminant-like hits =",
    len(discriminant_candidates)
)

for ell, a, label, expr in discriminant_candidates:

    print(
        f"  ell={ell}, a={a}, {label}:"
    )

    print(
        f"    {expr}"
    )

print()


# ============================================================================
# [8] NORMALIZATION BY N
# ============================================================================

print("[8] NORMALIZED N-ONLY SEARCH")
print("-" * 78)

normalized_hits = []

for ell, a, label, expr in discriminant_candidates:

    for denom_name, denom in [
        ("N", N),
        ("N^2", N**2),
        ("N*(N+1)", N*(N+1)),
        ("N*(N+2)", N*(N+2)),
    ]:

        candidate = sp.cancel(
            expr / denom
        )

        derivative = sp.diff(
            candidate,
            S
        )

        if zero(
            sp.together(
                derivative
            )
        ):

            normalized_hits.append(
                (
                    ell,
                    a,
                    label,
                    denom_name,
                    fact(candidate)
                )
            )

print(
    "  normalized pure-N hits =",
    len(normalized_hits)
)

for item in normalized_hits[:30]:

    ell, a, label, denom_name, expr = item

    print(
        f"  ell={ell}, a={a}, {label}, /{denom_name}"
    )

    print(
        f"    {expr}"
    )

print()


# ============================================================================
# [9] DIRECT RECOVERY OF S AND N FROM A PAIR
# ============================================================================

print("[9] TWO-POINT RECOVERY CERTIFICATE")
print("-" * 78)

# For ell=2, the exact shifted kernel is linear:
#
# H_(1,2)(x) = 6N - S^2 + S*x
#
# Hence:
#
#   H(1)-H(-1) = 2S
#
# and after recovering S:
#
#   H(1)+H(-1) = 12N - 2S^2
#
# allowing N to be recovered.

H2x = H_polys[2]

Hp = sp.expand(
    H2x.subs(
        x,
        1
    )
)

Hm = sp.expand(
    H2x.subs(
        x,
        -1
    )
)

S_pair = sp.cancel(
    (Hp - Hm) / 2
)

N_pair = sp.cancel(
    (
        (Hp + Hm)
        + 2*S_pair**2
    ) / 12
)

print(
    "  H_(1,2)(1)  =",
    fact(Hp)
)

print(
    "  H_(1,2)(-1) =",
    fact(Hm)
)

print(
    "  recovered S =",
    fact(S_pair)
)

print(
    "  S difference =",
    fact(
        S_pair - S
    )
)

print(
    "  recovered N =",
    fact(N_pair)
)

print(
    "  N difference =",
    fact(
        N_pair - N
    )
)

pair_recovery_pass = (
    zero(S_pair - S)
    and zero(N_pair - N)
)

print(
    "  TWO-POINT RECOVERY PASS =",
    pair_recovery_pass
)

print()


# ============================================================================
# [10] CONNECTION TO SHIFTED PRODUCT
# ============================================================================

print("[10] SHIFTED-PRODUCT CONNECTION")
print("-" * 78)

M1 = sp.expand(
    N + S + 1
)

M1_from_pair = sp.expand(
    N_pair + S_pair + 1
)

print(
    "  true M1 =",
    M1
)

print(
    "  pair-recovered M1 =",
    fact(M1_from_pair)
)

print(
    "  difference =",
    fact(
        M1_from_pair - M1
    )
)

print(
    "  PASS =",
    zero(
        M1_from_pair - M1
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

    # H_(1,2)(x) = 6N - S^2 + S*x
    hp = 6*nn - ss**2 + ss
    hm = 6*nn - ss**2 - ss

    recovered_s = (
        hp - hm
    ) // 2

    recovered_n = (
        hp
        + hm
        + 2*recovered_s**2
    ) // 12

    passed = (
        recovered_s == ss
        and recovered_n == nn
    )

    if not passed:
        numeric_failures += 1

    print(
        f"  ({pp},{qq})"
    )

    print(
        f"    N={nn}"
    )

    print(
        f"    S={ss}"
    )

    print(
        f"    H(+1)={hp}"
    )

    print(
        f"    H(-1)={hm}"
    )

    print(
        f"    recovered S={recovered_s}"
    )

    print(
        f"    recovered N={recovered_n}"
    )

    print(
        f"    PASS={passed}"
    )

print()

print(
    "  NUMERICAL FAILURES =",
    numeric_failures
)

print()


# ============================================================================
# [12] SEARCH FOR A SHIFTED-PRODUCT FROM PAIRED DATA
# ============================================================================

print("[12] PAIRED DATA -> SHIFTED PRODUCT SEARCH")
print("-" * 78)

# For any a:
#
# M_a = N + a*S + a^2.
#
# Using paired H_(1,2)(+a), H_(1,2)(-a):
#
# S = [H(+a)-H(-a)]/(2a).
#
# We test whether M_a can be written compactly.

for a in [1, 2]:

    Hp = sp.expand(
        H2x.subs(
            x,
            a
        )
    )

    Hm = sp.expand(
        H2x.subs(
            x,
            -a
        )
    )

    S_a = sp.cancel(
        (Hp - Hm) / (2*a)
    )

    M_a = sp.cancel(
        N + a*S_a + a**2
    )

    expected = (
        N + a*S + a**2
    )

    print(
        f"  a={a}:"
    )

    print(
        f"    recovered S =",
        fact(S_a)
    )

    print(
        f"    M_a =",
        fact(M_a)
    )

    print(
        f"    M_a difference =",
        fact(
            M_a - expected
        )
    )

    print(
        f"    PASS =",
        zero(
            M_a - expected
        )
    )

print()


# ============================================================================
# [13] DOES SYMMETRIC SHIFTING CREATE NEW N-ONLY INFORMATION?
# ============================================================================

print("[13] STRUCTURAL CONCLUSION")
print("-" * 78)

pure_pair_hits = len(
    discriminant_candidates
)

pure_quadratic_hits = len(
    quadratic_candidates
)

print(
    "  pure-N paired invariants =",
    pure_pair_hits
)

print(
    "  pure-N quadratic products =",
    pure_quadratic_hits
)

print()

print(
    "  Interpretation:"
)

print(
    "    The x -> -x symmetry can separate even and odd"
)

print(
    "    components of the shifted kernel."
)

print(
    "    For k=1, the odd component directly contains S,"
)

print(
    "    while the even component contains the remaining"
)

print(
    "    combination of N and S^2."
)

print(
    "    This gives an exact two-evaluation coordinate system"
)

print(
    "    for (N,S), but does not by itself provide an N-only"
)

print(
    "    source-free construction."
)

print()


# ============================================================================
# [14] FINAL STATUS
# ============================================================================

print("[14] EXPERIMENT STATUS")
print("-" * 78)

all_core = (
    all(
        zero(
            sp.diff(
                expr,
                S
            )
        )
        for expr in [
            fact(
                H_polys[2].subs(
                    x,
                    1
                )
                - (
                    6*N
                    - S**2
                    + S
                )
            ),
            fact(
                H_polys[2].subs(
                    x,
                    -1
                )
                - (
                    6*N
                    - S**2
                    - S
                )
            ),
        ]
    )
)

print(
    "  shifted-kernel construction =",
    all_core
)

print(
    "  paired S recovery =",
    pair_recovery_pass
)

print(
    "  numerical recovery =",
    numeric_failures == 0
)

print(
    "  pure-N paired invariants found =",
    pure_pair_hits
)

print(
    "  pure-N quadratic products found =",
    pure_quadratic_hits
)

print()

print(
    "MAIN RESULT:"
)

print(
    "  The x -> -x shifted-kernel symmetry gives:"
)

print(
    "      odd part  -> S"
)

print(
    "      even part -> 6N-S^2"
)

print(
    "  and therefore two evaluations H(+a), H(-a)"
)

print(
    "  recover both N and S exactly."
)

print()

print(
    "BUT:"
)

print(
    "  this still assumes the shifted kernel itself"
)

print(
    "  can be generated independently."
)

print()

print(
    "NEXT BRIDGE:"
)

print(
    "  Determine whether the paired shifted evaluations"
)

print(
    "  can be constructed from the existing N-only"
)

print(
    "  KAPPA/homogeneous-layer objects."
)

print()

print("=" * 78)
print("EXPERIMENT 492 FINISHED")
print("=" * 78)
