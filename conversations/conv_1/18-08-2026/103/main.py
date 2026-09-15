#!/usr/bin/env python3

import math
import sympy as sp

k, p = sp.symbols("k p")


# =============================================================================
# EXPERIMENT 106
# EXACT PARITY-SECTOR INDEX BOUNDARY LAW AUDIT
# =============================================================================
#
# This is a corrected, self-contained version.
#
# Main established structure:
#
#     centered parity coefficient c_p(k)
#       = (K-k)_(s_p) Q_p(k)
#
# We test whether s_p itself has a simple exact law in p, and inspect the
# resulting quotient Q_p(k).
#
# Arithmetic is exact QQ.
# Floating point is forbidden.
# =============================================================================


# =============================================================================
# DATA FROM EXPERIMENT 104
# =============================================================================

A_even = {
    0: [-sp.Rational(12879), -sp.Rational(28241), -sp.Rational(26989),
        -sp.Rational(13611), -sp.Rational(17875, 6),
        -sp.Rational(116923, 1680), -sp.Rational(5, 144)],
    2: [sp.Rational(2797337, 1920), sp.Rational(8986567, 2880),
        sp.Rational(36298273, 13440), sp.Rational(11670379, 11520),
        sp.Rational(19954213, 161280), -sp.Rational(9389, 4032), 0],
    4: [-sp.Rational(2083937, 30720), -sp.Rational(83529, 640),
        -sp.Rational(7965025, 129024), sp.Rational(2225141, 46080),
        sp.Rational(710501, 215040), 0, 0],
    6: [sp.Rational(85591, 61440), sp.Rational(83651, 46080),
        -sp.Rational(3174439, 2580480), 0, 0, 0, 0],
    8: [-sp.Rational(4913, 491520), 0, 0, 0, 0, 0, 0],
}

A_odd = {
    1: [sp.Rational(12143, 560), sp.Rational(2699231, 1344),
        sp.Rational(9120441, 2240), sp.Rational(21975383, 6720),
        sp.Rational(4460869, 5760), sp.Rational(42929, 1680), 0],
    3: [-sp.Rational(989, 11520), -sp.Rational(12024227, 46080),
        -sp.Rational(175956721, 322560), -sp.Rational(67903883, 161280),
        -sp.Rational(2590159, 53760), 0, 0],
    5: [-sp.Rational(517, 23040), sp.Rational(396119, 30720),
        sp.Rational(26625517, 1290240), -sp.Rational(234707, 129024),
        0, 0, 0],
    7: [sp.Rational(373, 1290240), -sp.Rational(1028053, 5160960),
        0, 0, 0, 0, 0],
}

B_even = {
    0: [sp.Rational(12980463, 1024), sp.Rational(6255583, 256),
        sp.Rational(87841139, 4480), sp.Rational(16998339, 2240),
        sp.Rational(1091983, 896), sp.Rational(1553, 240)],
    2: [-sp.Rational(19344659, 15360), -sp.Rational(26986999, 11520),
        -sp.Rational(2066529, 1120), -sp.Rational(573325, 576),
        sp.Rational(3312053, 40320), 0],
    4: [sp.Rational(129415, 3072), sp.Rational(267779, 3840),
        sp.Rational(224417, 4480), -sp.Rational(101119, 5040), 0, 0],
    6: [-sp.Rational(2267, 5120), -sp.Rational(6053, 11520),
        0, 0, 0, 0],
}

B_odd = {
    1: [-sp.Rational(584531, 35840), -sp.Rational(2908483, 1920),
        -sp.Rational(31233169, 13440), -sp.Rational(1446167, 1344),
        -sp.Rational(22259149, 40320), -sp.Rational(301, 240)],
    3: [-sp.Rational(59257, 46080), sp.Rational(186547, 1440),
        sp.Rational(367433, 1680), sp.Rational(126549, 448),
        -sp.Rational(162139, 40320), 0],
    5: [sp.Rational(4457, 46080), -sp.Rational(16819, 5760),
        -sp.Rational(5769, 896), 0, 0, 0],
    7: [-sp.Rational(421, 322560), 0, 0, 0, 0, 0],
}


CHANNELS = {
    "A-even": (A_even, 6),
    "A-odd": (A_odd, 6),
    "B-even": (B_even, 5),
    "B-odd": (B_odd, 5),
}


# =============================================================================
# EXACT HELPERS
# =============================================================================

def clean(expr):
    return sp.factor(sp.cancel(sp.expand(sp.sympify(expr))))


def is_zero(expr):
    return clean(expr) == 0


def poly(expr, var):
    return sp.Poly(
        sp.expand(sp.sympify(expr)),
        var,
        domain=sp.QQ,
    )


def degree(expr, var):
    expr = clean(expr)
    if expr == 0:
        return -1
    return int(poly(expr, var).degree())


def coefficient(expr, var, power):
    """
    Exact coefficient [var^power] expr.

    This was the missing helper that caused the previous NameError.
    """
    expr = clean(expr)

    if expr == 0:
        return sp.Rational(0)

    P = poly(expr, var)

    return sp.Rational(
        P.nth(int(power))
    )


def falling(x, n):
    out = sp.Integer(1)

    for r in range(int(n)):
        out *= x - r

    return clean(out)


def interpolate_values(values, var):
    points = [
        (sp.Integer(i), sp.Rational(v))
        for i, v in enumerate(values)
    ]

    return clean(
        sp.interpolate(points, var)
    )


def exact_division(numerator, denominator, var):
    """
    Exact polynomial division over QQ.
    Returns None when there is a nonzero remainder.
    """
    N = poly(numerator, var)
    D = poly(denominator, var)

    if D.is_zero:
        raise ZeroDivisionError("Polynomial divisor is zero.")

    q, r = sp.div(
        N,
        D,
        domain=sp.QQ,
    )

    if not r.is_zero:
        return None

    return clean(q.as_expr())


def primitive_integer_signature(expr, var):
    """
    Convert rational polynomial coefficients into primitive integer
    coefficients, highest power first.

    Handles constants and zero polynomials safely.
    """
    expr = clean(expr)

    if expr == 0:
        return []

    P = poly(expr, var)
    coeffs = [
        sp.Rational(c)
        for c in P.all_coeffs()
    ]

    if not coeffs:
        return []

    lcm_den = 1

    for c in coeffs:
        lcm_den = math.lcm(
            lcm_den,
            int(c.q),
        )

    ints = [
        int(c * lcm_den)
        for c in coeffs
    ]

    g = 0
    for value in ints:
        g = math.gcd(
            g,
            abs(value),
        )

    if g:
        ints = [
            value // g
            for value in ints
        ]

    for value in ints:
        if value != 0 and value < 0:
            ints = [-v for v in ints]
            break

    return ints


def support_profile(values):
    nz = [
        i for i, value in enumerate(values)
        if sp.Rational(value) != 0
    ]

    K = len(values) - 1

    if not nz:
        return {
            "k_min": None,
            "k_max": None,
            "tail": K + 1,
        }

    k_max = max(nz)

    return {
        "k_min": min(nz),
        "k_max": k_max,
        "tail": K - k_max,
    }


def exact_polynomial_law(records, max_degree):
    """
    Find an exact polynomial s(p) of degree <= max_degree
    through every supplied point.
    """
    if len(records) < max_degree + 1:
        return None

    for deg_target in range(max_degree + 1):

        if len(records) < deg_target + 1:
            continue

        initial = records[:deg_target + 1]

        P = clean(
            sp.interpolate(
                [
                    (
                        sp.Integer(pp),
                        sp.Integer(ss),
                    )
                    for pp, ss in initial
                ],
                p,
            )
        )

        if degree(P, p) > deg_target:
            continue

        ok = all(
            clean(
                P.subs(p, pp) - ss
            ) == 0
            for pp, ss in records
        )

        if ok:
            return P

    return None


# =============================================================================
# 1. SUPPORT / TERMINAL MULTIPLICITY
# =============================================================================

print("=" * 78)
print("EXPERIMENT 106 — EXACT PARITY-SECTOR INDEX BOUNDARY LAW AUDIT")
print("=" * 78)

print()
print("=" * 78)
print("1. EXACT SUPPORT / TERMINAL MULTIPLICITY")
print("=" * 78)

profiles = {}

for name, (table, K) in CHANNELS.items():

    print()
    print(name)

    records = []

    for power in sorted(table):

        info = support_profile(
            table[power]
        )

        records.append(
            (
                sp.Integer(power),
                sp.Integer(info["tail"]),
            )
        )

        print(
            f"  y^{power}: "
            f"k_max={info['k_max']} "
            f"K={K} "
            f"tail_s={info['tail']}"
        )

    profiles[name] = records


# =============================================================================
# 2. MULTIPLICITY LAWS
# =============================================================================

print()
print("=" * 78)
print("2. EXACT MULTIPLICITY LAWS IN CENTERED POWER")
print("=" * 78)

for name, records in profiles.items():

    print()
    print(name)

    found = False

    for deg_target in range(0, 5):

        P = exact_polynomial_law(
            records,
            deg_target,
        )

        if P is not None:

            print(
                f"  degree<={deg_target}: EXACT"
            )

            print(
                f"    s(p) = {P}"
            )

            found = True
            break

    if not found:
        print(
            "  no exact polynomial law through degree 4"
        )


# =============================================================================
# 3. DIRECT CENTERED-POWER FORM CHECKS
# =============================================================================

print()
print("=" * 78)
print("3. DIRECT SUPPORT-LAW CHECKS")
print("=" * 78)

for name, records in profiles.items():

    print()
    print(name)

    for power, s in records:

        print(
            f"  y^{power}: "
            f"s={s}"
        )


# =============================================================================
# 4. EXACT TERMINAL FACTOR EXTRACTION
# =============================================================================

print()
print("=" * 78)
print("4. EXACT TERMINAL FACTOR EXTRACTION")
print("=" * 78)

quotients = {}

for name, (table, K) in CHANNELS.items():

    print()
    print(name)

    quotients[name] = {}

    for power in sorted(table):

        values = table[power]

        info = support_profile(values)
        s = info["tail"]

        P = interpolate_values(
            values,
            k,
        )

        divisor = falling(
            K - k,
            s,
        )

        Q = exact_division(
            P,
            divisor,
            k,
        )

        ok = Q is not None

        if ok:
            quotients[name][power] = Q

        print(
            f"  y^{power}: "
            f"s={s} "
            f"exact={ok}"
        )

        if ok:
            print(
                f"    divisor = {divisor}"
            )
            print(
                f"    Q(k) = {Q}"
            )
            print(
                f"    degree(Q) = {degree(Q, k)}"
            )


# =============================================================================
# 5. QUOTIENT-DEGREE LAW
# =============================================================================

print()
print("=" * 78)
print("5. QUOTIENT DEGREE VS CENTERED POWER")
print("=" * 78)

for name in CHANNELS:

    print()
    print(name)

    records = []

    for power in sorted(quotients[name]):

        Q = quotients[name][power]

        dq = degree(
            Q,
            k,
        )

        records.append(
            (
                sp.Integer(power),
                sp.Integer(dq),
            )
        )

        print(
            f"  y^{power}: "
            f"degree(Q)={dq}"
        )

    found = False

    for deg_target in range(0, 5):

        P = exact_polynomial_law(
            records,
            deg_target,
        )

        if P is not None:

            print(
                f"  exact degree law "
                f"degree<={deg_target}: "
                f"deg(Q)={P}"
            )

            found = True
            break

    if not found:
        print(
            "  no polynomial degree law through degree 4"
        )


# =============================================================================
# 6. QUOTIENT SIGNATURES
# =============================================================================

print()
print("=" * 78)
print("6. TERMINAL-QUOTIENT PRIMITIVE SIGNATURES")
print("=" * 78)

for name in CHANNELS:

    print()
    print(name)

    for power in sorted(quotients[name]):

        Q = quotients[name][power]

        sig = primitive_integer_signature(
            Q,
            k,
        )

        print(
            f"  y^{power}: "
            f"{sig}"
        )


# =============================================================================
# 7. SECOND-STAGE LOWER-END FACTOR SEARCH
# =============================================================================

print()
print("=" * 78)
print("7. SECOND-STAGE INDEX FACTOR SEARCH")
print("=" * 78)

for name in CHANNELS:

    print()
    print(name)

    for power in sorted(quotients[name]):

        Q = quotients[name][power]

        found = []

        qdeg = degree(
            Q,
            k,
        )

        # k_(r)
        for r in range(
            1,
            min(7, qdeg + 1),
        ):

            divisor = falling(
                k,
                r,
            )

            if exact_division(
                Q,
                divisor,
                k,
            ) is not None:
                found.append(
                    f"k_({r})"
                )

        # (k-a)_r
        for a in range(1, 4):

            for r in range(
                1,
                min(6, qdeg + 1),
            ):

                divisor = falling(
                    k - a,
                    r,
                )

                if exact_division(
                    Q,
                    divisor,
                    k,
                ) is not None:
                    found.append(
                        f"(k-{a})_({r})"
                    )

        print(
            f"  y^{power}: "
            f"{found}"
        )


# =============================================================================
# 8. CROSS-PARITY SUPPORT LAW
# =============================================================================

print()
print("=" * 78)
print("8. CROSS-PARITY SUPPORT-LAW COMPARISON")
print("=" * 78)

for left, right in [
    ("A-even", "B-even"),
    ("A-odd", "B-odd"),
]:

    L = dict(profiles[left])
    R = dict(profiles[right])

    common = sorted(
        set(L) & set(R)
    )

    print()
    print(
        f"{left} vs {right}"
    )

    for power in common:

        equal = (
            L[power] == R[power]
        )

        print(
            f"  y^{power}: "
            f"{left}={L[power]} "
            f"{right}={R[power]} "
            f"equal={equal}"
        )


# =============================================================================
# 9. CROSS-PARITY QUOTIENT COMPARISON
# =============================================================================

print()
print("=" * 78)
print("9. CROSS-PARITY TERMINAL-QUOTIENT COMPARISON")
print("=" * 78)

for left, right in [
    ("A-even", "B-even"),
    ("A-odd", "B-odd"),
]:

    L = quotients[left]
    R = quotients[right]

    common = sorted(
        set(L) & set(R)
    )

    print()
    print(
        f"{left} vs {right}"
    )

    for power in common:

        qL = L[power]
        qR = R[power]

        same = (
            clean(qL - qR) == 0
        )

        ratio = None

        if not same:

            for kk in range(0, 7):

                a = clean(
                    qL.subs(k, kk)
                )

                b = clean(
                    qR.subs(k, kk)
                )

                if b != 0:

                    ratio = clean(
                        a / b
                    )
                    break

        print(
            f"  y^{power}: "
            f"same_exact={same} "
            f"first_nonzero_ratio={ratio}"
        )


# =============================================================================
# 10. EXACT POINTWISE RECONSTRUCTION
# =============================================================================

print()
print("=" * 78)
print("10. EXACT RECONSTRUCTION")
print("=" * 78)

reconstruction = True

for name, (table, K) in CHANNELS.items():

    print()
    print(name)

    for power in sorted(table):

        values = table[power]

        info = support_profile(values)
        s = info["tail"]

        divisor = falling(
            K - k,
            s,
        )

        P = interpolate_values(
            values,
            k,
        )

        Q = quotients[name].get(
            power
        )

        if Q is None:

            reconstruction = False

            print(
                f"  y^{power}: exact=False"
            )

            continue

        ok = all(
            clean(
                (
                    divisor * Q
                ).subs(k, idx)
                - sp.Rational(value)
            ) == 0
            for idx, value
            in enumerate(values)
        )

        reconstruction &= ok

        print(
            f"  y^{power}: exact={ok}"
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
  Experiment 105 established exact terminal index factors for the
  centered parity coefficients:

      c_p(k) = (K-k)_(s_p) Q_p(k).

  This experiment isolates the boundary mechanism.

  The first object is the support boundary:

      s_p = K - max{k : c_p(k) != 0}.

  The second object is the residual magnitude:

      Q_p(k) = c_p(k) / (K-k)_(s_p).

  The important distinction is:

      SUPPORT LAW
          explains where coefficients become exactly zero;

      QUOTIENT LAW
          explains the nonzero coefficients inside that support.

  A common exact law for s_p across the four parity sectors would
  be strong evidence that the terminal zeros are produced by a
  genuine indexed boundary mechanism.

  Failure of a second-stage k-factor search means that the terminal
  factor may already be the natural complete index-side factor.

  All calculations here are exact over QQ.

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

print(
    f"  terminal_factor_extraction = "
    f"{all(len(quotients[name]) == len(table)
           for name, (table, _) in CHANNELS.items())}"
)

print(
    f"  reconstruction = {reconstruction}"
)

print(
    f"  failures = {0 if reconstruction else 1}"
)

print(
    f"  ALL BASIC CHECKS PASS = {reconstruction}"
)

print()
print("EXPERIMENT 106 COMPLETE")