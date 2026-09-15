#!/usr/bin/env python3

import sympy as sp

print("EXPERIMENT 496 START")
print("=" * 78)
print("CROSS-SHIFT TRANSFER: x=1 KAPPA -> x=-1 KERNEL")
print("=" * 78)
print()

# ============================================================================
# SYMBOLS
# ============================================================================

p, q, x = sp.symbols("p q x")
N, S = sp.symbols("N S")


def zero(expr):
    return sp.expand(expr) == 0


def fact(expr):
    return sp.factor(sp.expand(expr))


def rat(expr):
    return sp.cancel(sp.together(expr))


def sym_reduce(expr):
    """
    Exact symmetric reduction p,q -> S=p+q, N=pq.
    """
    expr = sp.together(expr)

    num, den = sp.fraction(expr)

    num = sp.expand(num)
    den = sp.expand(den)

    num_red, num_rem, num_map = sp.symmetrize(
        num,
        [p, q],
        formal=True,
    )

    den_red, den_rem, den_map = sp.symmetrize(
        den,
        [p, q],
        formal=True,
    )

    if num_rem != 0:
        raise ValueError(
            f"Numerator is not symmetric: {num_rem}"
        )

    if den_rem != 0:
        raise ValueError(
            f"Denominator is not symmetric: {den_rem}"
        )

    s1_num = num_map[0][0]
    s2_num = num_map[1][0]

    s1_den = den_map[0][0]
    s2_den = den_map[1][0]

    num_red = num_red.subs({
        s1_num: S,
        s2_num: N,
    })

    den_red = den_red.subs({
        s1_den: S,
        s2_den: N,
    })

    return sp.cancel(
        sp.together(
            num_red / den_red
        )
    )


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
        # Compute the exact polynomial quotient first.
        X = sp.symbols("X")
        F = F_shifted(k, ell, X)
        return sp.expand(
            sp.cancel(F / X).subs(X, 0)
        )

    return sp.expand(
        sp.cancel(
            F_shifted(k, ell, xv) / xv
        )
    )


# ============================================================================
# [1] ORDINARY KAPPA VALUES
# ============================================================================

print("[1] ORIGINAL x=1 KAPPA SEQUENCE")
print("-" * 78)

F = {}

for ell in range(2, 11):

    raw = H_shifted(
        1,
        ell,
        1,
    )

    F[ell] = sym_reduce(
        raw
    )

    print(
        f"  F_{ell} = {fact(F[ell])}"
    )

print()


# ============================================================================
# [2] NEGATIVE-SHIFT KERNEL
# ============================================================================

print("[2] x=-1 SHIFTED KERNEL")
print("-" * 78)

Hminus = {}

for ell in range(2, 11):

    raw = H_shifted(
        1,
        ell,
        -1,
    )

    Hminus[ell] = sym_reduce(
        raw
    )

    print(
        f"  H_-1({ell}) = {fact(Hminus[ell])}"
    )

print()


# ============================================================================
# [3] DIRECT SHIFT-DIFFERENCE IDENTITIES
# ============================================================================

print("[3] x=+1 / x=-1 DIFFERENCE")
print("-" * 78)

for ell in range(2, 9):

    plus = F[ell]
    minus = Hminus[ell]

    difference = fact(
        sp.expand(
            plus - minus
        )
    )

    print(
        f"  ell={ell}:"
    )

    print(
        f"    H(+1)-H(-1) = {difference}"
    )

print()


# ============================================================================
# [4] SEARCH FOR PURE-N SHIFT DIFFERENCES
# ============================================================================

print("[4] PURE-N SHIFT-DIFFERENCE SEARCH")
print("-" * 78)

pure_n_hits = []

for ell in range(2, 11):

    difference = sp.expand(
        F[ell] - Hminus[ell]
    )

    if sp.diff(
        difference,
        S
    ) == 0:

        pure_n_hits.append(
            (
                ell,
                fact(difference)
            )
        )

        print(
            f"  ell={ell}: {fact(difference)}"
        )

print(
    "  PURE-N HITS =",
    len(pure_n_hits)
)

print()


# ============================================================================
# [5] EXACT TARGETS 12*a*N
# ============================================================================

print("[5] NORMALIZED SHIFT DIFFERENCES")
print("-" * 78)

for ell, expr in pure_n_hits:

    quot_N = sp.cancel(
        expr / N
    )

    print(
        f"  ell={ell}: difference/N = {fact(quot_N)}"
    )

print()


# ============================================================================
# [6] TARGET CROSS-SHIFT TRANSFER
# ============================================================================

print("[6] CROSS-SHIFT TRANSFER TARGET")
print("-" * 78)

print(
    "The key target is:"
)

print(
    "  H_(1,3)(+1) - H_(1,3)(-1) = 12*N"
)

target_difference = sp.expand(
    F[3] - Hminus[3]
)

print(
    "  computed difference =",
    fact(target_difference)
)

print(
    "  target difference   = 12*N"
)

print(
    "  PASS =",
    zero(
        target_difference - 12*N
    )
)

print()


# ============================================================================
# [7] LINEAR TRANSFER SEARCH
# ============================================================================

print("[7] LINEAR x=1 -> x=-1 TRANSFER SEARCH")
print("-" * 78)

print(
    "Search for:"
)

print(
    "  H_-1(ell) = sum_i P_i(N) * F_i"
)

print(
    "where P_i(N) has bounded degree."
)

print()

# Search spans F_2 ... F_m with coefficient degree 0..4.
transfer_hits = []

for max_ell in range(3, 11):

    source = [
        F[i]
        for i in range(
            2,
            max_ell + 1
        )
    ]

    for coeff_degree in range(0, 5):

        coeff_symbols = []

        for i in range(
            2,
            max_ell + 1
        ):

            coeff_symbols.extend(
                [
                    sp.Symbol(
                        f"c_{i}_{d}"
                    )
                    for d in range(
                        coeff_degree + 1
                    )
                ]
            )

        coeff_map = {}
        cursor = 0

        coeffs = {}

        for i in range(
            2,
            max_ell + 1
        ):

            polynomial = 0

            for d in range(
                coeff_degree + 1
            ):

                symbol = coeff_symbols[
                    cursor
                ]

                cursor += 1

                polynomial += (
                    symbol * N**d
                )

            coeffs[i] = polynomial

        for target_ell in range(
            3,
            max_ell + 1
        ):

            target = Hminus[
                target_ell
            ]

            expression = sp.expand(
                sum(
                    coeffs[i] * F[i]
                    for i in range(
                        2,
                        max_ell + 1
                    )
                )
                - target
            )

            poly = sp.Poly(
                expression,
                S,
                N,
            )

            equations = []

            for monom, coeff in poly.terms():

                equations.append(
                    coeff
                )

            equations = [
                eq
                for eq in equations
                if eq != 0
            ]

            if not equations:
                continue

            matrix, rhs = sp.linear_eq_to_matrix(
                equations,
                coeff_symbols
            )

            solution = sp.linsolve(
                (
                    matrix,
                    rhs
                ),
                coeff_symbols
            )

            if solution == sp.EmptySet:
                continue

            # Determine whether a solution exists.
            solution_list = list(
                solution
            )

            if not solution_list:
                continue

            point = solution_list[0]

            # Reject families containing free symbols:
            # we are interested in concrete transfer coefficients.
            free = set()

            for value in point:

                free.update(
                    value.free_symbols
                )

            free.difference_update(
                set(
                    coeff_symbols
                )
            )

            if free:
                continue

            candidate = sp.expand(
                sum(
                    point[
                        coeff_symbols.index(
                            sp.Symbol(
                                f"c_{i}_{d}"
                            )
                        )
                    ]
                    * N**d
                    for i in range(
                        2,
                        max_ell + 1
                    )
                    for d in range(
                        coeff_degree + 1
                    )
                )
            )

            transfer_hits.append(
                (
                    target_ell,
                    max_ell,
                    coeff_degree,
                    point
                )
            )

            print(
                f"  HIT: target ell={target_ell}, "
                f"source through={max_ell}, "
                f"degree={coeff_degree}"
            )

print(
    "  TRANSFER HITS =",
    len(transfer_hits)
)

print()


# ============================================================================
# [8] SPECIAL TWO-VALUE SEARCH FOR H_-1(3)
# ============================================================================

print("[8] MINIMAL TRANSFER SEARCH FOR H_-1(3)")
print("-" * 78)

target = Hminus[3]

minimal_hits = []

for degree in range(0, 7):

    c2 = sp.symbols(
        f"a0:{degree+1}"
    )

    c3 = sp.symbols(
        f"b0:{degree+1}"
    )

    P2 = sum(
        c2[d] * N**d
        for d in range(
            degree + 1
        )
    )

    P3 = sum(
        c3[d] * N**d
        for d in range(
            degree + 1
        )
    )

    difference = sp.expand(
        P2 * F[2]
        + P3 * F[3]
        - target
    )

    poly = sp.Poly(
        difference,
        S,
        N,
    )

    equations = [
        coeff
        for monom, coeff in poly.terms()
        if coeff != 0
    ]

    variables = list(c2) + list(c3)

    matrix, rhs = sp.linear_eq_to_matrix(
        equations,
        variables,
    )

    sol = sp.linsolve(
        (
            matrix,
            rhs
        ),
        variables,
    )

    if sol == sp.EmptySet:
        continue

    solutions = list(sol)

    if not solutions:
        continue

    point = solutions[0]

    free_symbols = set()

    for value in point:
        free_symbols.update(
            value.free_symbols
        )

    free_symbols -= set(
        variables
    )

    if free_symbols:
        continue

    candidate = sp.expand(
        sum(
            point[d] * N**d
            for d in range(
                degree + 1
            )
        )
    )

    minimal_hits.append(
        (
            degree,
            point
        )
    )

    print(
        f"  degree <= {degree}: EXACT HIT"
    )

print(
    "  MINIMAL TRANSFER HITS =",
    len(minimal_hits)
)

print()


# ============================================================================
# [9] DIRECT RECONSTRUCTION OF THE PURE-N DIFFERENCE
# ============================================================================

print("[9] RECONSTRUCT PURE-N DIFFERENCE FROM x=1 KERNEL")
print("-" * 78)

print(
    "For ell=3:"
)

print(
    "  H(+1)-H(-1) = 12*N"
)

print(
    "Since H(+1)=F3,"
)

print(
    "  H(-1) = F3 - 12*N"
)

candidate_hminus3 = sp.expand(
    F[3] - 12*N
)

print(
    "  F3 - 12*N =",
    fact(candidate_hminus3)
)

print(
    "  target H_-1(3) =",
    fact(Hminus[3])
)

print(
    "  PASS =",
    zero(
        candidate_hminus3
        - Hminus[3]
    )
)

print()


# ============================================================================
# [10] CAN THE RESULT BE RELATED TO F2?
# ============================================================================

print("[10] RELATION TO F2")
print("-" * 78)

relation = sp.cancel(
    Hminus[3] / F[2]
)

print(
    "  H_-1(3) / F2 =",
    fact(relation)
)

print()


# ============================================================================
# [11] x=+1 / x=-1 RELATIONS THROUGH ell=2,3
# ============================================================================

print("[11] TWO-SHIFT SYSTEM")
print("-" * 78)

Hp2 = F[2]
Hm2 = Hminus[2]

Hp3 = F[3]
Hm3 = Hminus[3]

diff2 = sp.expand(
    Hp2 - Hm2
)

diff3 = sp.expand(
    Hp3 - Hm3
)

print(
    "  H(+1,2)-H(-1,2) =",
    fact(diff2)
)

print(
    "  H(+1,3)-H(-1,3) =",
    fact(diff3)
)

print()

print(
    "  diff2 / 2 =",
    fact(
        sp.cancel(
            diff2 / 2
        )
    )
)

print(
    "  diff3 / 12 =",
    fact(
        sp.cancel(
            diff3 / 12
        )
    )
)

print()


# ============================================================================
# [12] NUMERICAL AUDIT
# ============================================================================

print("[12] NUMERICAL AUDIT")
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

    hp2 = int(
        H_shifted(
            1,
            2,
            1
        ).subs(
            {
                p: pp,
                q: qq,
            }
        )
    )

    hm2 = int(
        H_shifted(
            1,
            2,
            -1
        ).subs(
            {
                p: pp,
                q: qq,
            }
        )
    )

    hp3 = int(
        H_shifted(
            1,
            3,
            1
        ).subs(
            {
                p: pp,
                q: qq,
            }
        )
    )

    hm3 = int(
        H_shifted(
            1,
            3,
            -1
        ).subs(
            {
                p: pp,
                q: qq,
            }
        )
    )

    recovered_S = (
        hp2 - hm2
    ) // 2

    recovered_N = (
        hp3 - hm3
    ) // 12

    local_pass = (
        recovered_S == ss
        and recovered_N == nn
        and hp3 - hm3 == 12 * nn
    )

    if not local_pass:
        numeric_failures += 1

    print(
        f"  ({pp},{qq})"
    )

    print(
        f"    N = {nn}"
    )

    print(
        f"    S = {ss}"
    )

    print(
        f"    H2(+1) = {hp2}"
    )

    print(
        f"    H2(-1) = {hm2}"
    )

    print(
        f"    H3(+1) = {hp3}"
    )

    print(
        f"    H3(-1) = {hm3}"
    )

    print(
        f"    recovered S = {recovered_S}"
    )

    print(
        f"    recovered N = {recovered_N}"
    )

    print(
        f"    PASS = {local_pass}"
    )

print()

print(
    "  NUMERICAL FAILURES =",
    numeric_failures
)

print()


# ============================================================================
# [13] PROOF CERTIFICATES
# ============================================================================

print("[13] SYMBOLIC PROOF CERTIFICATES")
print("-" * 78)

cert_1 = zero(
    H_shifted(1, 2, 1)
    - H_shifted(1, 2, -1)
    - 2*S
)

cert_2 = zero(
    H_shifted(1, 3, 1)
    - H_shifted(1, 3, -1)
    - 12*N
)

cert_3 = zero(
    H_shifted(1, 3, -1)
    - (
        H_shifted(1, 3, 1)
        - 12*N
    )
)

print(
    "  ell=2 difference certificate =",
    cert_1
)

print(
    "  ell=3 difference certificate =",
    cert_2
)

print(
    "  H_-1(3)=F3-12N certificate =",
    cert_3
)

print()


# ============================================================================
# [14] INFORMATION-MODEL CHECK
# ============================================================================

print("[14] INFORMATION-MODEL CHECK")
print("-" * 78)

print(
    "Established by this experiment:"
)

print(
    "  H(+1,2)-H(-1,2) = 2*S"
)

print(
    "  H(+1,3)-H(-1,3) = 12*N"
)

print(
    "Thus the pair of shifts x=+1 and x=-1"
)

print(
    "provides direct coordinates for (N,S)."
)

print()

print(
    "The critical unresolved question is whether"
)

print(
    "H(-1,ell) can be generated from the ordinary"
)

print(
    "x=1 KAPPA sequence using N-only coefficients."
)

print()


# ============================================================================
# [15] FINAL STATUS
# ============================================================================

print("[15] EXPERIMENT STATUS")
print("-" * 78)

print(
    "  pure-N shift differences found =",
    len(pure_n_hits)
)

print(
    "  exact ell=2 S bridge =",
    cert_1
)

print(
    "  exact ell=3 N bridge =",
    cert_2
)

print(
    "  H_-1(3) transfer identity =",
    cert_3
)

print(
    "  numerical audit =",
    numeric_failures == 0
)

print(
    "  linear transfer hits =",
    len(transfer_hits)
)

print(
    "  minimal H_-1(3) transfer hits =",
    len(minimal_hits)
)

print()

print(
    "MAIN TARGET"
)

print(
    "----------------------------------------------------------------------------"
)

print(
    "The experiment confirms an exact cross-shift coordinate system:"
)

print(
    "  x=+1 and x=-1"
)

print(
    "      -> S and N"
)

print()

print(
    "In particular:"
)

print(
    "  H_(1,3)(-1) = F_3 - 12*N."
)

print()

print(
    "Therefore a source-free construction of H_(1,3)(-1),"
)

print(
    "or an N-only method producing the same quantity,"
)

print(
    "would immediately close the remaining information gap."
)

print()

print("=" * 78)
print("EXPERIMENT 496 FINISHED")
print("=" * 78)
