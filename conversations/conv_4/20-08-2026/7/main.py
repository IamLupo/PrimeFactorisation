#!/usr/bin/env python3

from __future__ import annotations

import sympy as sp

print("EXPERIMENT 484 START")
print("=" * 78)
print("ORIGINAL KERNEL -> SHIFTED-PRODUCT / N-ONLY INVARIANT SEARCH")
print("=" * 78)
print()

# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")

# ============================================================================
# ORIGINAL KAPPA KERNEL, k=1
# ============================================================================

def F(ell: int):
    return sp.expand(
        p * (1 + q) ** ell
        + q * (1 + p) ** ell
        - p ** ell * (1 + q)
        - q ** ell * (1 + p)
    )


def symmetric_reduce(expr):
    expr = sp.expand(expr)

    reduced, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise AssertionError(
            f"Symmetrization remainder: {remainder}"
        )

    result = reduced

    for formal_var, concrete in mapping:
        if sp.expand(concrete - (p + q)) == 0:
            result = result.subs(formal_var, S)
        elif sp.expand(concrete - p * q) == 0:
            result = result.subs(formal_var, N)
        else:
            raise AssertionError(
                f"Unexpected symmetric mapping: {mapping}"
            )

    result = sp.expand(result)

    check = sp.expand(
        result.subs(
            {
                S: p + q,
                N: p * q,
            }
        ) - expr
    )

    if check != 0:
        raise AssertionError(
            "Symmetric reconstruction failed."
        )

    return result


# ============================================================================
# BUILD F_l
# ============================================================================

print("[1] ORIGINAL KAPPA VALUES")
print("-" * 78)

Fvals = {}

for ell in range(2, 11):
    Fvals[ell] = sp.factor(
        symmetric_reduce(F(ell))
    )
    print(
        f"  F_{ell} = {Fvals[ell]}"
    )

print()


# ============================================================================
# TARGET SHIFTED PRODUCTS
# ============================================================================

print("[2] TARGET SHIFTED PRODUCTS")
print("-" * 78)

print(
    "For x=1:"
)
print(
    "  M_1 = (p+1)(q+1) = N + S + 1"
)

print(
    "For general x:"
)
print(
    "  M_x = N + S*x + x^2"
)

print()

M1 = sp.expand(
    N + S + 1
)

print(
    f"  M_1 = {M1}"
)

print()


# ============================================================================
# [3] KNOWN TWO-VALUE BRIDGE
# ============================================================================

print("[3] KNOWN TWO-VALUE BRIDGE")
print("-" * 78)

ratio_32 = sp.factor(
    sp.cancel(
        Fvals[3] / Fvals[2]
    )
)

print(
    f"  F_3/F_2 = {ratio_32}"
)

print(
    f"  F_3/F_2 - 1 = "
    f"{sp.factor(sp.cancel(ratio_32 - 1))}"
)

print()


# ============================================================================
# [4] CAN F3/F2 GIVE M1?
# ============================================================================

print("[4] SHIFTED PRODUCT FROM F3/F2")
print("-" * 78)

candidate_M1 = sp.factor(
    sp.cancel(
        N + Fvals[3] / Fvals[2]
    )
)

target_M1 = M1

difference_M1 = sp.factor(
    sp.cancel(
        candidate_M1 - target_M1
    )
)

print(
    "Candidate:"
)
print(
    "  N + F_3/F_2"
)

print(
    f"  = {candidate_M1}"
)

print()

print(
    "Target:"
)
print(
    f"  M_1 = {target_M1}"
)

print()

print(
    f"  Difference = {difference_M1}"
)

print(
    "  PASS =",
    difference_M1 == 0
)

print()


# ============================================================================
# [5] SEARCH SIMPLE RATIOS
# ============================================================================

print("[5] SIMPLE F_l RATIO SEARCH")
print("-" * 78)

ratio_hits = []

for a in range(2, 9):
    for b in range(a + 1, 10):

        ratio = sp.factor(
            sp.cancel(
                Fvals[b] / Fvals[a]
            )
        )

        targets = {
            "S": S,
            "S+1": S + 1,
            "N": N,
            "N+S+1": N + S + 1,
            "N+S": N + S,
            "2N+S": 2 * N + S,
        }

        for name, target in targets.items():

            if sp.cancel(
                ratio - target
            ) == 0:

                ratio_hits.append(
                    (a, b, name, ratio)
                )

for hit in ratio_hits:
    print(
        f"  F_{hit[1]}/F_{hit[0]} = "
        f"{hit[2]} -> {hit[3]}"
    )

if not ratio_hits:
    print(
        "  No direct simple ratio hits beyond the known F3/F2=S+1."
    )

print()


# ============================================================================
# [6] SEARCH LINEAR COMBINATIONS
# ============================================================================

print("[6] LINEAR COMBINATION SEARCH")
print("-" * 78)

targets = {
    "N": N,
    "S": S,
    "S+1": S + 1,
    "N+S+1": N + S + 1,
    "N+S": N + S,
    "2N+S": 2*N + S,
}

linear_hits = []

# Test very small integer coefficient pairs.
for a in range(-3, 4):
    for b in range(-3, 4):

        if a == 0 and b == 0:
            continue

        for i in range(2, 7):
            for j in range(i + 1, 8):

                expr = sp.expand(
                    a * Fvals[i]
                    + b * Fvals[j]
                )

                for name, target in targets.items():

                    # Only accept exact scalar multiples
                    # with small integer scalar.
                    for c in range(-10, 11):

                        if c == 0:
                            continue

                        if sp.expand(
                            expr - c * target
                        ) == 0:

                            linear_hits.append(
                                (
                                    i,
                                    j,
                                    a,
                                    b,
                                    c,
                                    name,
                                )
                            )

for hit in linear_hits:
    print(
        "  "
        f"a={hit[2]}, b={hit[3]}: "
        f"a*F_{hit[0]} + b*F_{hit[1]} "
        f"= {hit[4]}*{hit[5]}"
    )

print(
    f"  LINEAR HITS = {len(linear_hits)}"
)

print()


# ============================================================================
# [7] QUADRATIC COMBINATION SEARCH
# ============================================================================

print("[7] QUADRATIC INVARIANT SEARCH")
print("-" * 78)

print(
    "Search simple determinant-like combinations:"
)

quadratic_candidates = {
    "F2*F4-F3^2":
        sp.expand(
            Fvals[2] * Fvals[4]
            - Fvals[3] ** 2
        ),

    "F2*F5-F3*F4":
        sp.expand(
            Fvals[2] * Fvals[5]
            - Fvals[3] * Fvals[4]
        ),

    "F2*F6-F3*F5":
        sp.expand(
            Fvals[2] * Fvals[6]
            - Fvals[3] * Fvals[5]
        ),

    "F3*F5-F4^2":
        sp.expand(
            Fvals[3] * Fvals[5]
            - Fvals[4] ** 2
        ),

    "F3*F6-F4*F5":
        sp.expand(
            Fvals[3] * Fvals[6]
            - Fvals[4] * Fvals[5]
        ),
}

quadratic_targets = {
    "N": N,
    "S": S,
    "N+S+1": N + S + 1,
    "N*(N+S+1)": N * (N + S + 1),
    "(N+S+1)^2": (N + S + 1) ** 2,
}

for name, expr in quadratic_candidates.items():

    print(
        f"  {name}:"
    )

    for target_name, target in quadratic_targets.items():

        # Test exact proportionality by rational simplification.
        if expr == 0:
            ratio = sp.Integer(0)
        else:
            ratio = sp.factor(
                sp.cancel(
                    expr / target
                )
            )

        if ratio.is_number:
            print(
                f"    -> proportional to {target_name}: "
                f"ratio={ratio}"
            )

print()


# ============================================================================
# [8] HANKEL / MOMENT STRUCTURE
# ============================================================================

print("[8] HANKEL / MOMENT DETERMINANT SEARCH")
print("-" * 78)

print(
    "Because F_l is a finite exponential sequence,"
)
print(
    "Hankel determinants may expose elementary symmetric"
)
print(
    "functions of the four exponential bases."
)
print()

H2 = sp.factor(
    sp.expand(
        Fvals[2] * Fvals[4]
        - Fvals[3] ** 2
    )
)

H3 = sp.factor(
    sp.expand(
        sp.Matrix(
            [
                [Fvals[2], Fvals[3], Fvals[4]],
                [Fvals[3], Fvals[4], Fvals[5]],
                [Fvals[4], Fvals[5], Fvals[6]],
            ]
        ).det()
    )
)

print(
    f"  H_2 = {H2}"
)

print(
    f"  H_3 = {H3}"
)

print()


# ============================================================================
# [9] FACTOR HANKEL DETERMINANTS
# ============================================================================

print("[9] HANKEL FACTORIZATION")
print("-" * 78)

print(
    "Factored H_2:"
)

print(
    f"  {sp.factor(H2)}"
)

print()

print(
    "Factored H_3:"
)

print(
    f"  {sp.factor(H3)}"
)

print()


# ============================================================================
# [10] NORMALIZED HANKEL RATIOS
# ============================================================================

print("[10] NORMALIZED HANKEL RATIO SEARCH")
print("-" * 78)

H2_ratios = {
    "H2/F2^2":
        sp.factor(
            sp.cancel(
                H2 / Fvals[2] ** 2
            )
        ),

    "H2/F2/F3":
        sp.factor(
            sp.cancel(
                H2 / (Fvals[2] * Fvals[3])
            )
        ),

    "H2/F3^2":
        sp.factor(
            sp.cancel(
                H2 / Fvals[3] ** 2
            )
        ),
}

for name, expr in H2_ratios.items():
    print(
        f"  {name} = {expr}"
    )

print()


# ============================================================================
# [11] CAN HANKEL OBJECTS EXPOSE M1?
# ============================================================================

print("[11] SHIFTED-PRODUCT TARGET TEST")
print("-" * 78)

H_targets = {
    "N": N,
    "S": S,
    "M1": N + S + 1,
    "N*M1": N * (N + S + 1),
}

H2_target_hits = []

for name, expr in H2_ratios.items():

    for target_name, target in H_targets.items():

        difference = sp.cancel(
            expr - target
        )

        if difference == 0:
            H2_target_hits.append(
                (name, target_name)
            )
            print(
                f"  HIT: {name} = {target_name}"
            )

if not H2_target_hits:
    print(
        "  No direct normalized Hankel target hit."
    )

print()


# ============================================================================
# [12] EXPONENTIAL BASE STRUCTURE
# ============================================================================

print("[12] FOUR-BASE STRUCTURE")
print("-" * 78)

print(
    "The k=1 sequence is:"
)
print(
    "  F_l = p(q+1)^l + q(p+1)^l"
)
print(
    "        - (q+1)p^l - (p+1)q^l"
)
print()

print(
    "Bases:"
)
print(
    "  q+1, p+1, p, q"
)
print()

print(
    "Their elementary symmetric data are:"
)

bases = [
    p,
    q,
    p + 1,
    q + 1,
]

e1 = sp.expand(
    sum(bases)
)

e2 = sp.expand(
    sum(
        bases[i] * bases[j]
        for i in range(4)
        for j in range(i + 1, 4)
    )
)

e3 = sp.expand(
    sum(
        bases[i] * bases[j] * bases[k]
        for i in range(4)
        for j in range(i + 1, 4)
        for k in range(j + 1, 4)
    )
)

e4 = sp.expand(
    sp.prod(bases)
)

print(
    f"  e1 = {symmetric_reduce(e1)}"
)

print(
    f"  e2 = {symmetric_reduce(e2)}"
)

print(
    f"  e3 = {symmetric_reduce(e3)}"
)

print(
    f"  e4 = {symmetric_reduce(e4)}"
)

print()


# ============================================================================
# [13] NEW OBSERVATION TARGET: M1
# ============================================================================

print("[13] NEW BRIDGE TARGET")
print("-" * 78)

print(
    "The shifted product is:"
)

print(
    "  M1 = (p+1)(q+1) = N+S+1"
)

print()

print(
    "Experiment 483 gives:"
)

print(
    "  S+1 = F3/F2"
)

print()

print(
    "Therefore:"
)

print(
    "  M1 = N + F3/F2"
)

print(
    "This is an exact identity whenever F2 != 0."
)

candidate = sp.factor(
    sp.cancel(
        N + Fvals[3] / Fvals[2]
    )
)

target = sp.expand(
    N + S + 1
)

candidate_pass = (
    sp.cancel(
        candidate - target
    ) == 0
)

print(
    f"  candidate = {candidate}"
)

print(
    f"  target    = {target}"
)

print(
    f"  PASS = {candidate_pass}"
)

print()


# ============================================================================
# [14] INFORMATION-MODEL ANALYSIS
# ============================================================================

print("[14] INFORMATION-MODEL ANALYSIS")
print("-" * 78)

print(
    "The exact chain established so far is:"
)

print()
print(
    "  F2,F3"
)
print(
    "    -> S+1 = F3/F2"
)
print(
    "    -> S = F3/F2 - 1"
)
print(
    "    -> M1 = N + F3/F2"
)
print(
    "    -> z^2-S*z+N"
)
print(
    "    -> p,q"
)
print()

print(
    "The unresolved step remains:"
)

print(
    "  N"
)
print(
    "   -> F2,F3"
)

print()

print(
    "A successful N-only construction must therefore generate"
)
print(
    "at least one genuinely new quantity equivalent to F2,F3,"
)
print(
    "or directly generate M1=N+S+1."
)

print()


# ============================================================================
# [15] PROOF CERTIFICATES
# ============================================================================

print("[15] PROOF CERTIFICATES")
print("-" * 78)

cert_1 = sp.expand(
    Fvals[3] - (S + 1) * Fvals[2]
) == 0

cert_2 = sp.cancel(
    N + Fvals[3] / Fvals[2]
    - (N + S + 1)
) == 0

print(
    f"  F3-(S+1)F2 = 0: {cert_1}"
)

print(
    f"  N+F3/F2-(N+S+1) = 0: {cert_2}"
)

print()


# ============================================================================
# [16] STATUS
# ============================================================================

overall = (
    cert_1
    and cert_2
    and candidate_pass
)

print("[16] EXPERIMENT STATUS")
print("-" * 78)

print(
    f"  F3=(S+1)F2 exact: {cert_1}"
)

print(
    f"  M1=N+F3/F2 exact: {cert_2}"
)

print(
    f"  Shifted-product identity recovered: "
    f"{candidate_pass}"
)

print(
    f"  OVERALL EXACT AUDIT = {overall}"
)

print()

print(
    "MAIN RESULT:"
)

print(
    "  The original KAPPA kernel gives:"
)

print(
    "      F3/F2 = S+1"
)

print(
    "  and therefore:"
)

print(
    "      N + F3/F2 = (p+1)(q+1)"
)

print()

print(
    "This reconnects the modern KAPPA kernel directly"
)

print(
    "to the shifted-product identity explored in 2023."
)

print()

print(
    "NEXT TARGET:"
)

print(
    "  Determine whether (p+1)(q+1)=N+S+1 can itself"
)

print(
    "  be generated from the existing N-only constructions"
)

print(
    "  without first generating F2 or F3."
)

print()

print("=" * 78)
print("EXPERIMENT 484 FINISHED")
print("=" * 78)
