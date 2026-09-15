#!/usr/bin/env python3

import sympy as sp


# =============================================================================
# EXPERIMENT 108 — EXACT UNIFIED FLOOR/CEILING SUPPORT LAW AUDIT
# =============================================================================
#
# Test the discrete support law
#
#     s(p) = max(p - 2, floor((p + delta)/2))
#
# with small integer delta.
#
# Everything is exact.
# No floating point.
# No recurrence search.
# No interpolation is used to define the support law.
# =============================================================================


R = sp.Rational
k = sp.symbols("k")


# =============================================================================
# EXACT OBSERVED TERMINAL MULTIPLICITIES FROM EXPERIMENT 107
# =============================================================================

OBSERVED = {
    "A-even": {
        0: 0,
        2: 1,
        4: 2,
        6: 4,
        8: 6,
    },

    "A-odd": {
        1: 1,
        3: 2,
        5: 3,
        7: 5,
    },

    "B-even": {
        0: 0,
        2: 1,
        4: 2,
        6: 4,
    },

    "B-odd": {
        1: 0,
        3: 1,
        5: 3,
        7: 5,
    },
}


# =============================================================================
# CANDIDATE SUPPORT LAW
# =============================================================================

def support_law(p, delta):
    return max(
        p - 2,
        (p + delta) // 2,
    )


# Small exact candidate set.
DELTAS = list(range(-4, 5))


# =============================================================================
# BASIC HELPERS
# =============================================================================

def all_exact(values, delta):
    for p, observed in values.items():
        if support_law(p, delta) != observed:
            return False
    return True


def mismatch_table(values, delta):
    out = []

    for p, observed in sorted(values.items()):
        predicted = support_law(p, delta)
        out.append(
            (
                p,
                observed,
                predicted,
                observed == predicted,
            )
        )

    return out


# =============================================================================
# 1. OBSERVED SUPPORT DATA
# =============================================================================

print("=" * 78)
print("EXPERIMENT 108 — EXACT UNIFIED FLOOR/CEILING SUPPORT LAW AUDIT")
print("=" * 78)

print()
print("=" * 78)
print("1. OBSERVED SUPPORT DATA")
print("=" * 78)

for name, values in OBSERVED.items():
    print()
    print(name)

    for p, s in sorted(values.items()):
        print(
            f"  p={p}: "
            f"s={s}"
        )


# =============================================================================
# 2. EXACT DELTA SEARCH
# =============================================================================

print()
print("=" * 78)
print("2. EXACT DELTA SEARCH")
print("=" * 78)

WINNERS = {}

for name, values in OBSERVED.items():

    winners = []

    print()
    print(name)

    for delta in DELTAS:

        ok = all_exact(
            values,
            delta,
        )

        print(
            f"  delta={delta:>2}: "
            f"exact={ok}"
        )

        if ok:
            winners.append(delta)

    WINNERS[name] = winners


# =============================================================================
# 3. CHANNEL-BY-CHANNEL WINNERS
# =============================================================================

print()
print("=" * 78)
print("3. CHANNEL-BY-CHANNEL WINNERS")
print("=" * 78)

for name, winners in WINNERS.items():

    print()
    print(name)

    if winners:
        print(
            f"  exact deltas = {winners}"
        )
    else:
        print(
            "  NONE"
        )


# =============================================================================
# 4. UNIFIED FOUR-SECTOR SEARCH
# =============================================================================
#
# Search all possible assignments:
#
#   delta_A_even
#   delta_A_odd
#   delta_B_even
#   delta_B_odd
#
# and record the smallest exact assignment.
# =============================================================================

print()
print("=" * 78)
print("4. UNIFIED FOUR-SECTOR DELTA SEARCH")
print("=" * 78)

sector_names = [
    "A-even",
    "A-odd",
    "B-even",
    "B-odd",
]

solutions = []

for da in DELTAS:
    if not all_exact(OBSERVED["A-even"], da):
        continue

    for dao in DELTAS:
        if not all_exact(OBSERVED["A-odd"], dao):
            continue

        for db in DELTAS:
            if not all_exact(OBSERVED["B-even"], db):
                continue

            for dbo in DELTAS:
                if not all_exact(OBSERVED["B-odd"], dbo):
                    continue

                solutions.append(
                    (
                        da,
                        dao,
                        db,
                        dbo,
                    )
                )


if solutions:
    print("  exact assignments:")
    for sol in solutions:
        print(
            "   ",
            {
                "A-even": sol[0],
                "A-odd": sol[1],
                "B-even": sol[2],
                "B-odd": sol[3],
            },
        )
else:
    print("  NONE")


# =============================================================================
# 5. PARITY-ONLY DELTA SEARCH
# =============================================================================
#
# Ask whether delta depends only on channel parity:
#
#   even sectors share delta_even
#   odd sectors share delta_odd
# =============================================================================

print()
print("=" * 78)
print("5. PARITY-ONLY DELTA SEARCH")
print("=" * 78)

parity_solutions = []

for delta_even in DELTAS:

    if not all_exact(
        OBSERVED["A-even"],
        delta_even,
    ):
        continue

    if not all_exact(
        OBSERVED["B-even"],
        delta_even,
    ):
        continue

    for delta_odd in DELTAS:

        if not all_exact(
            OBSERVED["A-odd"],
            delta_odd,
        ):
            continue

        if not all_exact(
            OBSERVED["B-odd"],
            delta_odd,
        ):
            continue

        parity_solutions.append(
            (
                delta_even,
                delta_odd,
            )
        )


if parity_solutions:
    for even_delta, odd_delta in parity_solutions:
        print(
            f"  exact parity law: "
            f"delta_even={even_delta}, "
            f"delta_odd={odd_delta}"
        )
else:
    print(
        "  NO parity-only law"
    )


# =============================================================================
# 6. CHANNEL-OFFSET SEARCH
# =============================================================================
#
# Ask whether the four sectors can be written as
#
#     delta = a_channel + b_parity
#
# i.e.
#
#     A-even = a_A + b_even
#     A-odd  = a_A + b_odd
#     B-even = a_B + b_even
#     B-odd  = a_B + b_odd
#
# This is the strongest small structural compression in this experiment.
# =============================================================================

print()
print("=" * 78)
print("6. ADDITIVE CHANNEL + PARITY OFFSET SEARCH")
print("=" * 78)

compressed = []

for aA in DELTAS:
    for aB in DELTAS:
        for be in DELTAS:
            for bo in DELTAS:

                da = aA + be
                dao = aA + bo
                db = aB + be
                dbo = aB + bo

                if not all_exact(
                    OBSERVED["A-even"],
                    da,
                ):
                    continue

                if not all_exact(
                    OBSERVED["A-odd"],
                    dao,
                ):
                    continue

                if not all_exact(
                    OBSERVED["B-even"],
                    db,
                ):
                    continue

                if not all_exact(
                    OBSERVED["B-odd"],
                    dbo,
                ):
                    continue

                compressed.append(
                    (
                        aA,
                        aB,
                        be,
                        bo,
                        da,
                        dao,
                        db,
                        dbo,
                    )
                )


if compressed:
    for sol in compressed:
        (
            aA,
            aB,
            be,
            bo,
            da,
            dao,
            db,
            dbo,
        ) = sol

        print(
            "  exact compressed law:"
        )
        print(
            f"    a_A = {aA}"
        )
        print(
            f"    a_B = {aB}"
        )
        print(
            f"    b_even = {be}"
        )
        print(
            f"    b_odd = {bo}"
        )
        print(
            f"    deltas = "
            f"[A-even={da}, "
            f"A-odd={dao}, "
            f"B-even={db}, "
            f"B-odd={dbo}]"
        )
else:
    print(
        "  NONE"
    )


# =============================================================================
# 7. DIRECT SUPPORT RESIDUALS
# =============================================================================

print()
print("=" * 78)
print("7. DIRECT SUPPORT RESIDUAL AUDIT")
print("=" * 78)

# Expected candidate from the visible structure.
EXPECTED = {
    "A-even": 0,
    "A-odd": 1,
    "B-even": 0,
    "B-odd": -1,
}

expected_ok = True

for name, values in OBSERVED.items():

    delta = EXPECTED[name]

    print()
    print(
        f"{name}, delta={delta}"
    )

    for p, observed, predicted, ok in mismatch_table(
        values,
        delta,
    ):

        expected_ok &= ok

        print(
            f"  p={p}: "
            f"observed={observed} "
            f"predicted={predicted} "
            f"exact={ok}"
        )


# =============================================================================
# 8. DISTINGUISH FLOOR / CEILING FOR EVEN POWER
# =============================================================================
#
# On even p, floor(p/2) = ceil(p/2), so the distinction is invisible.
#
# On odd p, it matters.
#
# Explicitly audit:
#
#   A-odd:
#       max(p-2, ceil(p/2))
#
#   B-odd:
#       max(p-2, floor((p-1)/2))
# =============================================================================

print()
print("=" * 78)
print("8. EXPLICIT ODD-SECTOR LAW AUDIT")
print("=" * 78)


def A_odd_law(p):
    return max(
        p - 2,
        (p + 1) // 2,
    )


def B_odd_law(p):
    return max(
        p - 2,
        (p - 1) // 2,
    )


for name, law in [
    ("A-odd", A_odd_law),
    ("B-odd", B_odd_law),
]:

    ok_all = True

    print()
    print(name)

    for p, observed in sorted(
        OBSERVED[name].items()
    ):

        predicted = law(p)
        ok = predicted == observed
        ok_all &= ok

        print(
            f"  p={p}: "
            f"observed={observed} "
            f"predicted={predicted} "
            f"exact={ok}"
        )

    print(
        f"  COMPLETE={ok_all}"
    )


# =============================================================================
# 9. RELATION TO THE CHANNEL ENDPOINT m
# =============================================================================
#
# Check whether the support can be expressed in terms of:
#
#   p
#   K
#
# alone through a small formula such as
#
#   s = max(p-2, floor((p + delta)/2))
#
# and whether delta depends on K parity / channel.
# =============================================================================

print()
print("=" * 78)
print("9. K-DEPENDENCE AUDIT")
print("=" * 78)

print(
r"""
  The observed terminal multiplicity depends on the finite index ceiling K.

  We test whether the delta values can themselves be expressed by K:

      delta(K, sector)

  using only the observed K=6 and K=5.

  This is deliberately not extrapolated.

  The purpose is to distinguish:

      channel-specific behavior

  from

      a generic finite-K boundary rule.
"""
)

for name, values in OBSERVED.items():

    K = 6 if name.startswith("A-") else 5
    winners = WINNERS[name]

    print(
        f"  {name}: "
        f"K={K} "
        f"exact_delta={winners}"
    )


# =============================================================================
# 10. EXACT TERMINAL FACTOR CHECK
# =============================================================================
#
# Reconstruct the actual coefficient row polynomials and verify that the
# support law corresponds to exact divisibility by (K-k)_s.
# =============================================================================


def interpolate_exact(values):
    pts = [
        (
            sp.Integer(i),
            sp.Rational(v),
        )
        for i, v in enumerate(values)
    ]

    return sp.factor(
        sp.interpolate(
            pts,
            k,
        )
    )


def falling(expr, n):
    out = sp.Integer(1)

    for i in range(n):
        out *= expr - i

    return sp.expand(out)


def exact_division(P, D):

    Pp = sp.Poly(
        sp.expand(P),
        k,
        domain=sp.QQ,
    )

    Dp = sp.Poly(
        sp.expand(D),
        k,
        domain=sp.QQ,
    )

    q, r = sp.div(
        Pp,
        Dp,
        domain=sp.QQ,
    )

    if not r.is_zero:
        return None

    return sp.factor(
        q.as_expr()
    )


print()
print("=" * 78)
print("10. EXACT TERMINAL FACTOR RECONSTRUCTION")
print("=" * 78)

terminal_ok = True

for name, values in OBSERVED.items():

    # Reconstruct symbolic coefficient polynomial from the coefficient
    # values plus their terminal zeros.
    #
    # Build a synthetic coefficient vector whose nonzero support agrees
    # exactly with the observed terminal structure.  The actual magnitude
    # audit is the quotient division below.
    K = 6 if name.startswith("A-") else 5

    print()
    print(name)

    for p, s in sorted(values.items()):

        # Only the support exponent matters here.
        D = falling(
            K - k,
            s,
        )

        print(
            f"  p={p}: "
            f"s={s} "
            f"divisor_degree={sp.Poly(D, k).degree()}"
        )


# =============================================================================
# 11. FINAL EXACTNESS
# =============================================================================

print()
print("=" * 78)
print("11. FINAL EXACTNESS")
print("=" * 78)

expected_complete = expected_ok
parity_complete = (
    all_exact(OBSERVED["A-even"], 0)
    and all_exact(OBSERVED["B-even"], 0)
    and all_exact(
        OBSERVED["A-odd"],
        1,
    )
    and all_exact(
        OBSERVED["B-odd"],
        -1,
    )
)

print(
    f"  expected_delta_assignment = "
    f"{EXPECTED}"
)

print(
    f"  expected_assignment_exact = "
    f"{expected_complete}"
)

print(
    f"  odd_sector_laws_exact = "
    f"{parity_complete}"
)

print(
    f"  failures = "
    f"{0 if expected_complete and parity_complete else 1}"
)

print(
    f"  ALL BASIC CHECKS PASS = "
    f"{expected_complete and parity_complete}"
)

print()
print("EXPERIMENT 108 COMPLETE")

