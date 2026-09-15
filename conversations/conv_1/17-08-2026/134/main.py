# main.py

from math import comb


print("=" * 78)
print("EXPERIMENT 199")
print("DIRECT pq KERNEL -> QUOTIENT -> NEWTON-LAURENT AUDIT")
print("=" * 78)
print()
print("This experiment starts from the exact pq-level kernel")
print()
print(" F(k,l;p,q) = p^k(1+q)^l + q^k(1+p)^l")
print("              - p^l(1+q)^k - q^l(1+p)^k")
print()
print("and uses")
print()
print(" F = (p+q+1) Q(pq,p+q).")
print()
print("No guessed convolution and no affine candidate search.")
print()


# ============================================================================
# BASIC HELPERS
# ============================================================================

def C(n, r):
    if r < 0 or r > n:
        return 0
    return comb(n, r)


def add_poly(a, b):
    out = dict(a)

    for e, c in b.items():
        out[e] = out.get(e, 0) + c

        if out[e] == 0:
            del out[e]

    return out


def scale_poly(a, scalar):
    if scalar == 0:
        return {}

    return {
        e: scalar * c
        for e, c in a.items()
        if scalar * c != 0
    }


def mul_poly(a, b):
    out = {}

    for ea, ca in a.items():
        for eb, cb in b.items():
            e = ea + eb
            out[e] = out.get(e, 0) + ca * cb

    return {
        e: c
        for e, c in out.items()
        if c != 0
    }


def x_poly():
    # x = u + u^(-1)
    return {
        1: 1,
        -1: 1,
    }


X = x_poly()


# ============================================================================
# NEWTON / LAURENT BASIS
# ============================================================================

def Pi(j):
    """
    Laurent realization

        Pi_0 = 2
        Pi_j = u^j + u^(-j), j >= 1
    """

    j = abs(j)

    if j == 0:
        return {0: 2}

    return {
        j: 1,
        -j: 1,
    }


def coeff(poly, exponent):
    return poly.get(exponent, 0)


def positive_support(poly):
    return sorted(
        [e for e in poly.keys() if e >= 0],
        reverse=True,
    )


def format_laurent(poly):
    if not poly:
        return "0"

    parts = []

    for e in sorted(poly.keys(), reverse=True):
        c = poly[e]

        if e == 0:
            parts.append(str(c))

        elif e == 1:
            parts.append(f"{c}*u")

        elif e == -1:
            parts.append(f"{c}*u^-1")

        elif e > 1:
            parts.append(f"{c}*u^{e}")

        else:
            parts.append(f"{c}*u^({e})")

    return " + ".join(parts)


# ============================================================================
# EXACT t-LAYER OF F
# ============================================================================
#
# From the paper:
#
# [t^m]F(tu,tu^-1)
#
# = C(l,m-k) Pi_|m-2k|
#   - C(k,m-l) Pi_|m-2l|
#
# This is implemented directly, without symbolic expansion.


def F_layer(k, ell, m):

    part1 = scale_poly(
        Pi(abs(m - 2 * k)),
        C(ell, m - k),
    )

    part2 = scale_poly(
        Pi(abs(m - 2 * ell)),
        C(k, m - ell),
    )

    return add_poly(
        part1,
        scale_poly(part2, -1),
    )


# ============================================================================
# EXACT QUOTIENT LAYERS
# ============================================================================
#
# Since
#
#     F(t) = (1 + t*x) Q(t)
#
# and
#
#     F_m = Q_m + x Q_{m-1},
#
# we recover recursively
#
#     Q_m = F_m - x Q_{m-1}.
#
# This is the crucial exact reconstruction.


def quotient_layers(k, ell):

    max_m = k + ell

    F_layers = [
        F_layer(k, ell, m)
        for m in range(max_m + 1)
    ]

    Q_layers = []

    previous = {}

    for m in range(max_m + 1):

        correction = mul_poly(X, previous)

        current = add_poly(
            F_layers[m],
            scale_poly(correction, -1),
        )

        Q_layers.append(current)
        previous = current

    return F_layers, Q_layers


# ============================================================================
# RECONSTRUCTION CHECK
# ============================================================================

def reconstruction_check(k, ell):

    F_layers, Q_layers = quotient_layers(k, ell)

    failures = 0

    for m in range(k + ell + 1):

        lhs = F_layers[m]

        previous = {}
        if m > 0:
            previous = Q_layers[m - 1]

        rhs = add_poly(
            Q_layers[m],
            mul_poly(X, previous),
        )

        if lhs != rhs:
            failures += 1

    return failures


# ============================================================================
# NEWTON COEFFICIENT EXTRACTION
# ============================================================================
#
# If a Laurent polynomial is symmetric
#
#     a_0 + a_1(u+u^-1) + a_2(u^2+u^-2) + ...
#
# then the coefficient of Pi_j is simply the coefficient of u^j.
#
# Therefore we can extract Newton coefficients directly.


def newton_coefficient(poly, j):

    j = abs(j)

    if j == 0:
        return coeff(poly, 0) // 2

    return coeff(poly, j)


# ============================================================================
# KNOWN PAPER CASES
# ============================================================================

CASES = [
    (3, 7),
    (3, 9),
    (5, 11),
    (5, 19),
    (7, 15),
    (9, 21),
    (11, 23),
]


# ============================================================================
# 1. EXACT KERNEL SANITY CHECK
# ============================================================================

print("=" * 78)
print("1. EXACT pq KERNEL SANITY CHECK")
print("=" * 78)

total_failures = 0

for k, ell in CASES:

    failures = reconstruction_check(k, ell)

    status = "PASS" if failures == 0 else "FAIL"

    print(
        f"k={k:2d} ell={ell:2d} "
        f"layer reconstruction failures={failures:3d} "
        f"{status}"
    )

    total_failures += failures

print()
print(
    f"total reconstruction failures = {total_failures}"
)


# ============================================================================
# 2. PRINT THE FIRST NONZERO F-LAYERS
# ============================================================================

print()
print("=" * 78)
print("2. EXACT F-LAYER STRUCTURE")
print("=" * 78)

for k, ell in CASES[:3]:

    print()
    print(f"k={k} ell={ell}")

    F_layers, _ = quotient_layers(k, ell)

    for m, layer in enumerate(F_layers):

        if not layer:
            continue

        print(
            f"  m={m:2d} : "
            f"{format_laurent(layer)}"
        )


# ============================================================================
# 3. QUOTIENT-LAYER STRUCTURE
# ============================================================================

print()
print("=" * 78)
print("3. QUOTIENT LAYERS")
print("=" * 78)

for k, ell in CASES[:4]:

    print()
    print(f"k={k} ell={ell}")

    _, Q_layers = quotient_layers(k, ell)

    for m, layer in enumerate(Q_layers):

        if not layer:
            continue

        print(
            f"  m={m:2d} : "
            f"{format_laurent(layer)}"
        )


# ============================================================================
# 4. NEWTON COEFFICIENT AUDIT
# ============================================================================
#
# For each Q_m inspect the coefficient of Pi_j.
#
# We do not yet claim that every j is the theorem's target coefficient.
# We simply expose the exact data.


print()
print("=" * 78)
print("4. NEWTON COEFFICIENT TABLE")
print("=" * 78)

for k, ell in CASES:

    print()
    print(f"k={k} ell={ell}")

    _, Q_layers = quotient_layers(k, ell)

    for m, layer in enumerate(Q_layers):

        if not layer:
            continue

        js = positive_support(layer)

        if not js:
            continue

        entries = []

        for j in js:
            value = newton_coefficient(layer, j)

            if value != 0:
                entries.append((j, value))

        if entries:

            print(
                f"  m={m:2d} : "
                + " ".join(
                    f"Pi_{j}={v}"
                    for j, v in entries
                )
            )


# ============================================================================
# 5. POST-BOUNDARY TARGET INDEX AUDIT
# ============================================================================
#
# The paper defines
#
#     j = ell - k - d - 1
#
# for d >= 1 and j >= 1.
#
# We inspect the corresponding exact coefficient in every Q-layer.


print()
print("=" * 78)
print("5. POST-BOUNDARY TARGET INDEX AUDIT")
print("=" * 78)

for k, ell in CASES:

    print()
    print(f"k={k} ell={ell}")

    _, Q_layers = quotient_layers(k, ell)

    found = 0

    for d in range(1, ell):

        j = ell - k - d - 1

        if j < 1:
            continue

        if j >= ell:
            continue

        values = []

        for m, layer in enumerate(Q_layers):

            a = newton_coefficient(layer, j)

            if a != 0:
                values.append((m, a))

        if not values:
            continue

        found += 1

        print(
            f"  d={d:2d} "
            f"j={j:2d} "
            f"coefficients={values}"
        )

    print(
        f"  target indices found = {found}"
    )


# ============================================================================
# 6. BOUNDARY-LAYER DIAGNOSTIC
# ============================================================================
#
# The paper distinguishes a top weighted layer and the next lower layer.
#
# We inspect the two highest nonzero Q layers directly.


print()
print("=" * 78)
print("6. TOP TWO QUOTIENT LAYERS")
print("=" * 78)

for k, ell in CASES:

    _, Q_layers = quotient_layers(k, ell)

    nonzero_m = [
        m
        for m, layer in enumerate(Q_layers)
        if layer
    ]

    print()
    print(f"k={k} ell={ell}")

    for m in nonzero_m[-2:]:

        layer = Q_layers[m]

        print(
            f"  m={m:2d} "
            f"support={positive_support(layer)}"
        )

        for j in positive_support(layer):

            a = newton_coefficient(layer, j)

            if a != 0:
                print(
                    f"      Pi_{j}: {a}"
                )


# ============================================================================
# 7. DIRECT COMPARISON WITH THE RECURRENCE CLAIMS
# ============================================================================
#
# Odd branch:
#
#     A_{j+2} + A_j = 0
#
# Even branch:
#
#     A_{j+4} + 2 A_{j+2} + A_j = 0
#
# We test these only inside the exact extracted layers.
#
# This is deliberately a recurrence audit, not a theorem claim.


print()
print("=" * 78)
print("7. RECURRENCE AUDIT")
print("=" * 78)


def recurrence_audit(layer, parity):

    values = {
        j: newton_coefficient(layer, j)
        for j in range(
            max(abs(min(layer.keys())), 1),
            max(abs(max(layer.keys())), 1) + 1
        )
    }

    failures = 0
    tests = 0

    j_values = sorted(
        j for j in values
        if j >= 1
    )

    if parity == "odd":
        for j in j_values:
            if j + 2 not in values:
                continue

            a = values[j]
            b = values[j + 2]

            tests += 1

            if a + b != 0:
                failures += 1

    else:
        for j in j_values:
            if (
                j + 2 not in values
                or j + 4 not in values
            ):
                continue

            a = values[j]
            b = values[j + 2]
            c = values[j + 4]

            tests += 1

            if c + 2 * b + a != 0:
                failures += 1

    return tests, failures


for k, ell in CASES:

    _, Q_layers = quotient_layers(k, ell)

    print()
    print(f"k={k} ell={ell}")

    for m, layer in enumerate(Q_layers):

        if not layer:
            continue

        odd_tests, odd_failures = recurrence_audit(
            layer,
            "odd",
        )

        even_tests, even_failures = recurrence_audit(
            layer,
            "even",
        )

        if odd_tests or even_tests:

            print(
                f"  m={m:2d} "
                f"odd={odd_failures}/{odd_tests} "
                f"even={even_failures}/{even_tests}"
            )


# ============================================================================
# 8. FINAL DIAGNOSTIC
# ============================================================================

print()
print("=" * 78)
print("FINAL DIAGNOSTIC")
print("=" * 78)
print()

print("This experiment has now crossed the critical boundary:")
print()
print("    p,q")
print("     |")
print("     v")
print("    F(k,l;p,q)")
print("     |")
print("     v")
print("    exact t-layers")
print("     |")
print("     v")
print("    Q(pq,p+q)")
print("     |")
print("     v")
print("    Newton coefficients")
print()

print(
    "No guessed pq -> (k,ell,s) map is used in the calculation."
)

print()
print(
    "NEXT MATHEMATICAL QUESTION:"
)
print(
    "Which exact Q-layer carries the theorem's post-boundary"
)
print(
    "leading N-coefficient, and does its Newton coefficient satisfy"
)
print(
    "the stated odd/even recurrence and boundary normalization?"
)

print()
print("=" * 78)
print("END EXPERIMENT 199")
print("=" * 78)

