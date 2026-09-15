#!/usr/bin/env python3

import sympy as sp

k = sp.symbols("k")


# =============================================================================
# EXPERIMENT 107 — EXACT PIECEWISE SUPPORT LAW / TERMINAL-FACTOR AUDIT
# =============================================================================
#
# Goal:
#
#   Experiment 106 found the terminal multiplicities
#
#       c_p(k) = (K-k)_(s_p) Q_p(k)
#
#   but used polynomial interpolation in p to describe s_p.
#
#   Here we test simple piecewise integer laws directly.
#
#   No interpolation is used for the support law.
#   All arithmetic is exact QQ.
#   Floating point is forbidden.
# =============================================================================


R = sp.Rational


# =============================================================================
# PARITY-COEFFICIENT TABLES
# =============================================================================

DATA = {

    "A-even": {
        "K": 6,
        "rows": {
            0: [R(-12879), R(-28241), R(-26989), R(-13611),
                R(-17875, 6), R(-116923, 1680), R(-5, 144)],

            2: [R(2797337, 1920), R(8986567, 2880),
                R(36298273, 13440), R(11670379, 11520),
                R(19954213, 161280), R(-9389, 4032), 0],

            4: [R(-2083937, 30720), R(-83529, 640),
                R(-7965025, 129024), R(2225141, 46080),
                R(710501, 215040), 0, 0],

            6: [R(85591, 61440), R(83651, 46080),
                R(-3174439, 2580480), 0, 0, 0, 0],

            8: [R(-4913, 491520), 0, 0, 0, 0, 0, 0],
        },
    },

    "A-odd": {
        "K": 6,
        "rows": {
            1: [R(12143, 560), R(2699231, 1344),
                R(9120441, 2240), R(21975383, 6720),
                R(4460869, 5760), R(42929, 1680), 0],

            3: [R(-989, 11520), R(-12024227, 46080),
                R(-175956721, 322560), R(-67903883, 161280),
                R(-2590159, 53760), 0, 0],

            5: [R(-517, 23040), R(396119, 30720),
                R(26625517, 1290240), R(-234707, 129024),
                0, 0, 0],

            7: [R(373, 1290240), R(-1028053, 5160960),
                0, 0, 0, 0, 0],
        },
    },

    "B-even": {
        "K": 5,
        "rows": {
            0: [R(12980463, 1024), R(6255583, 256),
                R(87841139, 4480), R(16998339, 2240),
                R(1091983, 896), R(1553, 240)],

            2: [R(-19344659, 15360), R(-26986999, 11520),
                R(-2066529, 1120), R(-573325, 576),
                R(3312053, 40320), 0],

            4: [R(129415, 3072), R(267779, 3840),
                R(224417, 4480), R(-101119, 5040), 0, 0],

            6: [R(-2267, 5120), R(-6053, 11520),
                0, 0, 0, 0],
        },
    },

    "B-odd": {
        "K": 5,
        "rows": {
            1: [R(-584531, 35840), R(-2908483, 1920),
                R(-31233169, 13440), R(-1446167, 1344),
                R(-22259149, 40320), R(-301, 240)],

            3: [R(-59257, 46080), R(186547, 1440),
                R(367433, 1680), R(126549, 448),
                R(-162139, 40320), 0],

            5: [R(4457, 46080), R(-16819, 5760),
                R(-5769, 896), 0, 0, 0],

            7: [R(-421, 322560), 0, 0, 0, 0, 0],
        },
    },
}


# =============================================================================
# BASIC EXACT HELPERS
# =============================================================================

def clean(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def poly(expr):
    return sp.Poly(
        sp.expand(expr),
        k,
        domain=sp.QQ,
    )


def degree(expr):
    expr = clean(expr)

    if expr == 0:
        return -1

    return int(poly(expr).degree())


def interpolate(values):
    pts = [
        (sp.Integer(i), sp.Rational(v))
        for i, v in enumerate(values)
    ]

    return clean(
        sp.interpolate(pts, k)
    )


def falling(expr, n):
    out = sp.Integer(1)

    for r in range(int(n)):
        out *= expr - r

    return clean(out)


def exact_quotient(P, D):
    Pp = poly(P)
    Dp = poly(D)

    if Dp.is_zero:
        raise ZeroDivisionError("zero polynomial divisor")

    q, r = sp.div(
        Pp,
        Dp,
        domain=sp.QQ,
    )

    if not r.is_zero:
        return None

    return clean(q.as_expr())


def actual_terminal_s(values):
    """
    s = K - largest index having a nonzero coefficient.
    """
    K = len(values) - 1

    nz = [
        i for i, value in enumerate(values)
        if sp.Rational(value) != 0
    ]

    if not nz:
        return K + 1

    return K - max(nz)


# =============================================================================
# CANDIDATE SUPPORT LAWS
# =============================================================================

def floor_half(p):
    return p // 2


def ceil_half(p):
    return (p + 1) // 2


def law_floor(p):
    return max(p // 2, p - 2)


def law_ceil(p):
    return max((p + 1) // 2, p - 2)


def law_linear(p):
    return max(0, p - 2)


def law_floor_only(p):
    return p // 2


def law_ceil_only(p):
    return (p + 1) // 2


CANDIDATES = {
    "max(floor(p/2), p-2)": law_floor,
    "max(ceil(p/2), p-2)": law_ceil,
    "max(0, p-2)": law_linear,
    "floor(p/2)": law_floor_only,
    "ceil(p/2)": law_ceil_only,
}


# =============================================================================
# 1. ACTUAL SUPPORT DATA
# =============================================================================

print("=" * 78)
print("EXPERIMENT 107 — EXACT PIECEWISE SUPPORT LAW / TERMINAL-FACTOR AUDIT")
print("=" * 78)

print()
print("=" * 78)
print("1. OBSERVED TERMINAL MULTIPLICITIES")
print("=" * 78)

actual = {}

for name, info in DATA.items():

    actual[name] = {}

    print()
    print(name)

    K = info["K"]

    for p, values in sorted(info["rows"].items()):

        s = actual_terminal_s(values)
        actual[name][p] = s

        print(
            f"  power={p}: "
            f"K={K} "
            f"s={s}"
        )


# =============================================================================
# 2. DIRECT PIECEWISE-LAW TEST
# =============================================================================

print()
print("=" * 78)
print("2. DIRECT PIECEWISE SUPPORT-LAW TEST")
print("=" * 78)

winning_laws = {}

for name in DATA:

    print()
    print(name)

    winners = []

    for law_name, law_fn in CANDIDATES.items():

        checks = []

        for p, observed in actual[name].items():

            predicted = law_fn(p)
            checks.append(
                predicted == observed
            )

        ok = all(checks)

        print(
            f"  {law_name}: "
            f"exact={ok}"
        )

        if ok:
            winners.append(law_name)

    winning_laws[name] = winners


# =============================================================================
# 3. BEST SIMPLE LAW BY CHANNEL
# =============================================================================

print()
print("=" * 78)
print("3. SIMPLEST EXACT SUPPORT LAW")
print("=" * 78)

for name, winners in winning_laws.items():

    print()
    print(name)

    if winners:
        print(
            f"  exact candidates = {winners}"
        )
    else:
        print(
            "  NONE"
        )


# =============================================================================
# 4. UNIFIED LAW AUDIT
# =============================================================================

print()
print("=" * 78)
print("4. UNIFIED LAW AUDIT")
print("=" * 78)

for law_name, law_fn in CANDIDATES.items():

    all_ok = True

    for name in DATA:
        for p, observed in actual[name].items():

            if law_fn(p) != observed:
                all_ok = False

    print(
        f"  {law_name}: "
        f"all_channels={all_ok}"
    )


# =============================================================================
# 5. EXACT TERMINAL FACTOR REMOVAL
# =============================================================================

print()
print("=" * 78)
print("5. EXACT TERMINAL FACTOR REMOVAL")
print("=" * 78)

quotients = {}

for name, info in DATA.items():

    K = info["K"]
    quotients[name] = {}

    print()
    print(name)

    for p, values in sorted(info["rows"].items()):

        s = actual[name][p]

        P = interpolate(values)

        # (K-k)_s
        D = falling(
            K - k,
            s,
        )

        Q = exact_quotient(
            P,
            D,
        )

        ok = Q is not None

        print(
            f"  power={p}: "
            f"s={s} "
            f"exact={ok}"
        )

        if not ok:
            continue

        quotients[name][p] = Q

        print(
            f"    divisor = {D}"
        )

        print(
            f"    quotient_degree = {degree(Q)}"
        )

        print(
            f"    Q(k) = {Q}"
        )


# =============================================================================
# 6. QUOTIENT DEGREE LAW
# =============================================================================

print()
print("=" * 78)
print("6. QUOTIENT DEGREE CHECK")
print("=" * 78)

for name, info in DATA.items():

    K = info["K"]

    print()
    print(name)

    for p, Q in sorted(
        quotients[name].items()
    ):

        observed = degree(Q)
        expected = K - actual[name][p]

        print(
            f"  power={p}: "
            f"observed={observed} "
            f"expected={expected} "
            f"exact={observed == expected}"
        )


# =============================================================================
# 7. TEST A COMMON FORMULA s(p) = max(h(p), p-2)
# =============================================================================
#
# h(p) is allowed to depend only on parity:
#
#   even p: a*p+b
#   odd  p: c*p+d
#
# with small half-integer-compatible rational coefficients.
#
# Since the data set is tiny, we enumerate small candidates exactly.
# =============================================================================

print()
print("=" * 78)
print("7. SMALL PIECEWISE-AFFINE BASE LAW SEARCH")
print("=" * 78)

# Base laws are represented as a*p + b, with a in {0,1/2,1}
# and b in {-1, -1/2, 0, 1/2, 1}.
SLOPES = [
    R(0),
    R(1, 2),
    R(1),
]

OFFSETS = [
    R(-1),
    R(-1, 2),
    R(0),
    R(1, 2),
    R(1),
]


def try_base_pair():

    results = []

    for ae in SLOPES:
        for be in OFFSETS:
            for ao in SLOPES:
                for bo in OFFSETS:

                    ok = True

                    for name in DATA:

                        for power, observed in actual[name].items():

                            if power % 2 == 0:
                                h = ae * power + be
                            else:
                                h = ao * power + bo

                            predicted = max(
                                h,
                                power - 2,
                            )

                            if predicted != observed:
                                ok = False
                                break

                        if not ok:
                            break

                    if ok:
                        results.append(
                            (ae, be, ao, bo)
                        )

    return results


pairs = try_base_pair()

if pairs:
    for ae, be, ao, bo in pairs:

        print(
            "  exact unified law:"
        )
        print(
            f"    even: h(p) = {ae}*p + {be}"
        )
        print(
            f"    odd : h(p) = {ao}*p + {bo}"
        )
else:
    print(
        "  NONE"
    )


# =============================================================================
# 8. SUPPORT LAW COMPARISON BETWEEN CHANNELS
# =============================================================================

print()
print("=" * 78)
print("8. CHANNEL SUPPORT COMPARISON")
print("=" * 78)

for pair in [
    ("A-even", "B-even"),
    ("A-odd", "B-odd"),
]:

    left, right = pair

    common = sorted(
        set(actual[left])
        & set(actual[right])
    )

    print()
    print(
        f"{left} vs {right}"
    )

    for p in common:

        a = actual[left][p]
        b = actual[right][p]

        print(
            f"  power={p}: "
            f"{left}={a} "
            f"{right}={b} "
            f"equal={a == b}"
        )


# =============================================================================
# 9. SECOND-STAGE FACTOR SEARCH ON QUOTIENTS
# =============================================================================

print()
print("=" * 78)
print("9. SECOND-STAGE EXACT k-FACTOR SEARCH")
print("=" * 78)

for name in quotients:

    print()
    print(name)

    for p, Q in sorted(
        quotients[name].items()
    ):

        found = []

        qdeg = degree(Q)

        # k_(r)
        for r in range(
            1,
            qdeg + 1,
        ):

            if exact_quotient(
                Q,
                falling(k, r),
            ) is not None:

                found.append(
                    f"k_({r})"
                )

        # (K-j)_r equivalent in the index k variable:
        # (K-k)_r
        K = DATA[name]["K"]

        for r in range(
            1,
            qdeg + 1,
        ):

            if exact_quotient(
                Q,
                falling(K - k, r),
            ) is not None:

                found.append(
                    f"({K}-k)_({r})"
                )

        print(
            f"  power={p}: {found}"
        )


# =============================================================================
# 10. EXACT RECONSTRUCTION
# =============================================================================

print()
print("=" * 78)
print("10. EXACT RECONSTRUCTION")
print("=" * 78)

all_reconstruct = True

for name, info in DATA.items():

    K = info["K"]

    print()
    print(name)

    for p, values in sorted(
        info["rows"].items()
    ):

        s = actual[name][p]
        Q = quotients[name][p]

        D = falling(
            K - k,
            s,
        )

        P_reconstructed = clean(
            D * Q
        )

        P_actual = interpolate(values)

        ok = clean(
            P_reconstructed - P_actual
        ) == 0

        all_reconstruct &= ok

        print(
            f"  power={p}: "
            f"exact={ok}"
        )


# =============================================================================
# 11. STRUCTURAL INTERPRETATION
# =============================================================================

print()
print("=" * 78)
print("11. STRUCTURAL INTERPRETATION")
print("=" * 78)

print(
r"""
  The previous experiment described terminal multiplicities with
  interpolated polynomials in the centered power p.

  That is not the right level of structure.

  The observed support is finite and integer-valued, so this experiment
  tests the more natural possibility

      s(p) = max( h_parity(p), p-2 )

  with a very small affine base law h_parity(p).

  The threshold term p-2 is motivated by the high-power rows, while
  the parity-dependent base captures the low-power boundary.

  A successful small piecewise-affine law is substantially stronger
  than a cubic interpolation because it explains the support pattern
  using a small discrete rule.

  After extracting the exact terminal factor

      (K-k)_(s(p)),

  the quotient Q_p(k) is examined independently.

  Thus the experiment separates:

      support geometry
          from
      coefficient magnitude.

  Everything is exact over QQ.

  No floating point.
  No recurrence search.
  No extrapolation.
"""
)


# =============================================================================
# 12. FINAL EXACTNESS
# =============================================================================

print()
print("=" * 78)
print("12. FINAL EXACTNESS")
print("=" * 78)

terminal_ok = True

for name, info in DATA.items():

    for p, values in info["rows"].items():

        s = actual[name][p]

        P = interpolate(values)

        D = falling(
            info["K"] - k,
            s,
        )

        Q = exact_quotient(
            P,
            D,
        )

        if Q is None:
            terminal_ok = False


print(
    f"  terminal_factor_extraction = {terminal_ok}"
)

print(
    f"  reconstruction = {all_reconstruct}"
)

print(
    f"  failures = {0 if terminal_ok and all_reconstruct else 1}"
)

print(
    f"  ALL BASIC CHECKS PASS = "
    f"{terminal_ok and all_reconstruct}"
)

print()
print("EXPERIMENT 107 COMPLETE")

