#!/usr/bin/env python3

import math
import sympy as sp

j, y, k = sp.symbols("j y k")


# =============================================================================
# EXPERIMENT 105
# EXACT PARITY-SECTOR INDEX FALLING-FACTORIAL AUDIT
# =============================================================================
#
# Established:
#
#   P_k(j)
#     = j_(k) (m-j)_(max(0,k-4)) R_k(j)
#
# and after centering
#
#   y = 2j-m
#
#   R_k(j) = E_k(y) + O_k(y).
#
# Experiment 104 showed strong terminal-zero patterns in the coefficients
# of E_k and O_k.
#
# Example:
#
#   A-even:
#       y^8 : nonzero only k=0
#       y^6 : nonzero only k=0,1,2
#       y^4 : nonzero only k=0,...,4
#
# This experiment asks whether those zeros are explained by an exact
# falling factorial in the INDEX k:
#
#     (K-k)_(s)
#
# where K is the largest observed channel index.
#
# Everything is exact over QQ.
# =============================================================================


# =============================================================================
# 1. ORIGINAL TRIANGULAR DATA
# =============================================================================

A = [
    [
        -2, -154, sp.Rational(-818), sp.Rational(-1360, 3),
        sp.Rational(8435, 24), sp.Rational(-5851, 120),
        sp.Rational(-13373, 720), sp.Rational(51773, 5040),
        sp.Rational(-4913, 1920),
    ],
    [
        -550, -5015, -4734, sp.Rational(7879, 3),
        sp.Rational(-5147, 120), sp.Rational(-210877, 720),
        sp.Rational(83651, 720), sp.Rational(-1028053, 40320),
    ],
    [
        -7125, -14567, 4635, sp.Rational(25508, 15),
        sp.Rational(-198919, 144), sp.Rational(427555, 1008),
        sp.Rational(-3174439, 40320),
    ],
    [
        -11900, -1711, sp.Rational(28949, 6),
        sp.Rational(-31711, 15), sp.Rational(404513, 840),
        sp.Rational(-234707, 4032),
    ],
    [
        sp.Rational(-17875, 6), sp.Rational(51337, 30),
        sp.Rational(-52447, 180), sp.Rational(-14333, 210),
        sp.Rational(710501, 13440),
    ],
    [
        sp.Rational(-1001, 12), sp.Rational(26687, 360),
        sp.Rational(-40921, 1260), sp.Rational(9389, 1008),
    ],
    [
        sp.Rational(-5, 72), sp.Rational(5, 72),
        sp.Rational(-5, 144),
    ],
]

B = [
    [
        25, 619, sp.Rational(3231, 2), sp.Rational(-33, 2),
        sp.Rational(-1675, 4), sp.Rational(3363, 20),
        sp.Rational(-9991, 360), sp.Rational(-421, 2520),
    ],
    [
        1750, 8624, sp.Rational(6829, 3),
        sp.Rational(-27341, 8), sp.Rational(10551, 10),
        sp.Rational(-16819, 180), sp.Rational(-6053, 180),
    ],
    [
        9690, 10234, sp.Rational(-57829, 8),
        sp.Rational(148151, 120), sp.Rational(1432, 5),
        sp.Rational(-5769, 28),
    ],
    [
        sp.Rational(22100, 3), sp.Rational(-19045, 12),
        sp.Rational(-5577, 4), sp.Rational(351271, 360),
        sp.Rational(-101119, 315),
    ],
    [
        sp.Rational(17875, 24), sp.Rational(-22061, 40),
        sp.Rational(132343, 720), sp.Rational(-162139, 5040),
    ],
    [
        sp.Rational(65, 12), sp.Rational(-313, 60),
        sp.Rational(301, 120),
    ],
]


# =============================================================================
# 2. SAFE EXACT HELPERS
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
    expr = clean(expr)
    if expr == 0:
        return sp.Integer(0)
    return sp.Rational(poly(expr, var).nth(power))


def falling(x, n):
    out = sp.Integer(1)
    for r in range(n):
        out *= x - r
    return clean(out)


def upper_falling(m, x, n):
    out = sp.Integer(1)
    for r in range(n):
        out *= m - x - r
    return clean(out)


def exact_division(num, den, var):
    num = poly(num, var)
    den = poly(den, var)

    q, r = sp.div(
        num,
        den,
        domain=sp.QQ,
    )

    if not r.is_zero:
        return None

    return clean(q.as_expr())


# =============================================================================
# 3. ORIGINAL CHANNEL RECONSTRUCTION
# =============================================================================

def reconstruct(row, k0):
    out = sp.Integer(0)

    for s, value in enumerate(row):
        r = k0 + s
        out += sp.sympify(value) * falling(j, r)

    return clean(out)


def residuals(channel, m):
    rows = []

    for k0, row in enumerate(channel):

        P = reconstruct(row, k0)

        s = max(0, k0 - 4)

        divisor = clean(
            falling(j, k0) *
            upper_falling(m, j, s)
        )

        R = exact_division(
            P,
            divisor,
            j,
        )

        if R is None:
            raise RuntimeError(
                f"Residual extraction failed at k={k0}"
            )

        rows.append(R)

    return rows


# =============================================================================
# 4. CENTERED PARITY
# =============================================================================

def centered(R, m):
    return clean(
        R.subs(
            j,
            (y + sp.Integer(m)) / 2,
        )
    )


def parity_parts(F):
    even = clean(
        (F + F.subs(y, -y)) / 2
    )

    odd = clean(
        (F - F.subs(y, -y)) / 2
    )

    return even, odd


# =============================================================================
# 5. EXACT PARITY COEFFICIENT TABLE
# =============================================================================

def parity_table(centered_rows, parity):
    """
    Dictionary:
        power -> [coefficient at k=0, ..., k=K]
    """
    maxdeg = max(
        [degree(R, y) for R in centered_rows] + [-1]
    )

    if maxdeg < 0:
        return {}

    if parity == "even":
        powers = list(range(0, maxdeg + 1, 2))
    else:
        powers = list(range(1, maxdeg + 1, 2))

    table = {}

    for power in powers:
        table[power] = [
            coefficient(
                R,
                y,
                power,
            )
            for R in centered_rows
        ]

    return table


# =============================================================================
# 6. DETECT TERMINAL ZERO BLOCK
# =============================================================================

def terminal_zero_count(values):
    """
    Number of consecutive zero values at the RIGHT end.
    """
    count = 0

    for value in reversed(values):
        if sp.Rational(value) == 0:
            count += 1
        else:
            break

    return count


# =============================================================================
# 7. CANDIDATE INDEX FALLING FACTOR
# =============================================================================

def candidate_index_factor(K, s):
    """
    (K-k)(K-k-1)...(K-k-s+1)
    """
    if s == 0:
        return sp.Integer(1)

    out = sp.Integer(1)

    for r in range(s):
        out *= K - k - r

    return clean(out)


# =============================================================================
# 8. EXACT INDEX DIVISIBILITY
# =============================================================================

def divide_index_factor(values, K, s):
    """
    Construct exact polynomial in k through the observed values,
    then test divisibility by the candidate index falling factorial.
    """
    if not values:
        return None, None

    P = clean(
        sp.interpolate(
            [
                (
                    sp.Integer(idx),
                    sp.Rational(value),
                )
                for idx, value in enumerate(values)
            ],
            k,
        )
    )

    divisor = candidate_index_factor(K, s)

    Q = exact_division(
        P,
        divisor,
        k,
    )

    return P, Q


# =============================================================================
# 9. PRIMITIVE INTEGER SIGNATURE
# =============================================================================

def primitive_signature(values):
    values = [
        sp.Rational(v)
        for v in values
    ]

    if not values:
        return []

    lcm_den = 1

    for value in values:
        lcm_den = math.lcm(
            lcm_den,
            int(value.q),
        )

    ints = [
        int(value * lcm_den)
        for value in values
    ]

    g = 0

    for value in ints:
        g = math.gcd(
            g,
            abs(value),
        )

    if g != 0:
        ints = [
            value // g
            for value in ints
        ]

    # Canonical sign.
    for value in ints:
        if value != 0:
            if value < 0:
                ints = [-z for z in ints]
            break

    return ints


# =============================================================================
# 10. MATRIX OF QUOTIENT COEFFICIENTS
# =============================================================================

def quotient_matrix(table, K):
    """
    For each parity power, remove the maximal terminal index factor
    and interpolate the quotient polynomial in k.

    Returns exact coefficient matrix.
    """
    quotients = []

    max_q_degree = -1

    for power, values in table.items():

        s = terminal_zero_count(values)

        P, Q = divide_index_factor(
            values,
            K,
            s,
        )

        if Q is None:
            # Keep a zero placeholder only if the data are identically zero.
            if all(
                sp.Rational(v) == 0
                for v in values
            ):
                Q = sp.Integer(0)
            else:
                raise RuntimeError(
                    f"Index falling-factor division failed for y^{power}"
                )

        quotients.append(
            (power, s, Q)
        )

        max_q_degree = max(
            max_q_degree,
            degree(Q, k),
        )

    if max_q_degree < 0:
        return sp.zeros(
            len(quotients),
            0,
        ), quotients

    data = []

    for power, s, Q in quotients:
        data.append([
            coefficient(
                Q,
                k,
                r,
            )
            for r in range(
                max_q_degree,
                -1,
                -1,
            )
        ])

    return sp.Matrix(data), quotients


# =============================================================================
# 11. RECONSTRUCTION FROM INDEX FACTORS
# =============================================================================

def verify_index_factor(table, K, power, s, Q):
    values = table[power]

    divisor = candidate_index_factor(
        K,
        s,
    )

    P = clean(
        Q * divisor
    )

    return all(
        clean(
            P.subs(
                k,
                idx,
            ) - sp.Rational(value)
        ) == 0
        for idx, value in enumerate(values)
    )


# =============================================================================
# 12. BUILD
# =============================================================================

A_R = residuals(A, 8)
B_R = residuals(B, 7)

A_centered = [
    centered(R, 8)
    for R in A_R
]

B_centered = [
    centered(R, 7)
    for R in B_R
]

A_even = [
    parity_parts(F)[0]
    for F in A_centered
]

A_odd = [
    parity_parts(F)[1]
    for F in A_centered
]

B_even = [
    parity_parts(F)[0]
    for F in B_centered
]

B_odd = [
    parity_parts(F)[1]
    for F in B_centered
]

tables = {
    "A-even": parity_table(A_even, "even"),
    "A-odd": parity_table(A_odd, "odd"),
    "B-even": parity_table(B_even, "even"),
    "B-odd": parity_table(B_odd, "odd"),
}

Ks = {
    "A-even": 6,
    "A-odd": 6,
    "B-even": 5,
    "B-odd": 5,
}


# =============================================================================
# OUTPUT
# =============================================================================

print("=" * 78)
print("EXPERIMENT 105 — EXACT PARITY-SECTOR INDEX FALLING-FACTORIAL AUDIT")
print("=" * 78)


# =============================================================================
# 1. TERMINAL ZERO PROFILE
# =============================================================================

print()
print("=" * 78)
print("1. TERMINAL-ZERO / INDEX-MULTIPLICITY PROFILE")
print("=" * 78)

for name, table in tables.items():

    K = Ks[name]

    print()
    print(name)

    for power, values in table.items():

        s = terminal_zero_count(values)

        print(
            f"  y^{power}: "
            f"values={values}"
        )

        print(
            f"      terminal_zero_count={s}"
        )

        print(
            f"      candidate_factor="
            f"{candidate_index_factor(K, s)}"
        )


# =============================================================================
# 2. EXACT INDEX FALLING DIVISIBILITY
# =============================================================================

print()
print("=" * 78)
print("2. EXACT INDEX FALLING DIVISIBILITY")
print("=" * 78)

all_divisions = True

quotient_data = {}

for name, table in tables.items():

    K = Ks[name]

    print()
    print(name)

    quotient_data[name] = []

    for power, values in table.items():

        s = terminal_zero_count(values)

        P, Q = divide_index_factor(
            values,
            K,
            s,
        )

        ok = Q is not None

        if not ok and all(
            sp.Rational(v) == 0
            for v in values
        ):
            Q = sp.Integer(0)
            ok = True

        all_divisions &= ok

        quotient_data[name].append(
            (power, s, P, Q)
        )

        print(
            f"  y^{power}: "
            f"exact={ok} "
            f"s={s}"
        )

        if ok:
            print(
                f"      Q(k)={Q}"
            )


# =============================================================================
# 3. QUOTIENT DEGREE PROFILE
# =============================================================================

print()
print("=" * 78)
print("3. INDEX-QUOTIENT DEGREE PROFILE")
print("=" * 78)

for name, rows in quotient_data.items():

    print()
    print(name)

    for power, s, P, Q in rows:

        print(
            f"  y^{power}: "
            f"original_degree={degree(P, k)} "
            f"factor_order={s} "
            f"quotient_degree={degree(Q, k)}"
        )


# =============================================================================
# 4. QUOTIENT PRIMITIVE SIGNATURES
# =============================================================================

print()
print("=" * 78)
print("4. QUOTIENT PRIMITIVE INTEGER SIGNATURES")
print("=" * 78)

for name, rows in quotient_data.items():

    print()
    print(name)

    for power, s, P, Q in rows:

        if Q == 0:
            sig = []

        else:
            d = degree(Q, k)

            coeffs = [
                coefficient(
                    Q,
                    k,
                    r,
                )
                for r in range(
                    d,
                    -1,
                    -1,
                )
            ]

            sig = primitive_signature(
                coeffs
            )

        print(
            f"  y^{power}: "
            f"{sig}"
        )


# =============================================================================
# 5. INDEX-QUOTIENT MATRIX RANK
# =============================================================================

print()
print("=" * 78)
print("5. INDEX-QUOTIENT MATRIX RANK")
print("=" * 78)

rank_data = {}

for name, table in tables.items():

    K = Ks[name]

    M, info = quotient_matrix(
        table,
        K,
    )

    rank_data[name] = (
        M,
        info,
    )

    print()
    print(
        f"  {name}: "
        f"shape={M.shape} "
        f"rank={M.rank()}"
    )


# =============================================================================
# 6. QUOTIENT MATRIX SEPARABILITY
# =============================================================================

print()
print("=" * 78)
print("6. QUOTIENT MATRIX SEPARABILITY")
print("=" * 78)

for name, (M, info) in rank_data.items():

    rank = M.rank()

    print()
    print(
        f"  {name}: "
        f"rank={rank}"
    )

    print(
        f"      rank-1={rank <= 1}"
    )


# =============================================================================
# 7. CROSS-PARITY COMPARISON
# =============================================================================

print()
print("=" * 78)
print("7. CROSS-PARITY INDEX-QUOTIENT COMPARISON")
print("=" * 78)

for left, right in [
    ("A-even", "B-even"),
    ("A-odd", "B-odd"),
]:

    print()
    print(
        f"{left} vs {right}"
    )

    M1, info1 = rank_data[left]
    M2, info2 = rank_data[right]

    min_rows = min(
        M1.rows,
        M2.rows,
    )

    min_cols = min(
        M1.cols,
        M2.cols,
    )

    for rr in range(min_rows):

        power1 = info1[rr][0]
        power2 = info2[rr][0]

        print(
            f"  row {rr}: "
            f"{left} y^{power1} vs "
            f"{right} y^{power2}"
        )

        ratios = []

        for cc in range(min_cols):

            a = M1[rr, cc]
            b = M2[rr, cc]

            if b == 0:
                continue

            if a == 0:
                ratios.append(
                    (cc, "0")
                )
            else:
                ratios.append(
                    (
                        cc,
                        sp.cancel(a / b),
                    )
                )

        print(
            f"      ratios={ratios}"
        )


# =============================================================================
# 8. EXACT RECONSTRUCTION AFTER INDEX FACTOR REMOVAL
# =============================================================================

print()
print("=" * 78)
print("8. EXACT RECONSTRUCTION AFTER INDEX FACTOR REMOVAL")
print("=" * 78)

reconstruction_ok = True

for name, rows in quotient_data.items():

    K = Ks[name]

    print()
    print(name)

    table = tables[name]

    for power, s, P, Q in rows:

        ok = verify_index_factor(
            table,
            K,
            power,
            s,
            Q,
        )

        reconstruction_ok &= ok

        print(
            f"  y^{power}: exact={ok}"
        )


# =============================================================================
# 9. SECOND-STAGE INDEX NORMALIZATION
# =============================================================================

print()
print("=" * 78)
print("9. SECOND-STAGE INDEX NORMALIZATION")
print("=" * 78)

print(
"""
  After removing the terminal factor

      (K-k)_s

  we inspect the remaining quotient Q(k).

  The next question is whether Q(k) itself contains a lower-end
  falling factorial k_(r), or a shifted factorial such as

      (k-a)_r.

  Only exact divisibility is accepted.
"""
)

for name, rows in quotient_data.items():

    print()
    print(name)

    for power, s, P, Q in rows:

        if Q == 0:
            print(
                f"  y^{power}: zero quotient"
            )
            continue

        candidates = []

        # Small candidate lower-end factors.
        for r in range(
            1,
            min(
                5,
                degree(Q, k) + 1,
            ),
        ):

            divisor = falling(
                k,
                r,
            )

            q = exact_division(
                Q,
                divisor,
                k,
            )

            if q is not None:
                candidates.append(
                    (
                        f"k_({r})",
                        q,
                    )
                )

        # Small shifted candidates.
        for shift in range(1, 4):

            for r in range(
                1,
                min(
                    4,
                    degree(Q, k) + 1,
                ),
            ):

                divisor = falling(
                    k - shift,
                    r,
                )

                q = exact_division(
                    Q,
                    divisor,
                    k,
                )

                if q is not None:
                    candidates.append(
                        (
                            f"(k-{shift})_({r})",
                            q,
                        )
                    )

        print(
            f"  y^{power}: "
            f"candidate_exact_factors={candidates}"
        )


# =============================================================================
# 10. FINAL INTERPRETATION
# =============================================================================

print()
print("=" * 78)
print("10. STRUCTURAL INTERPRETATION")
print("=" * 78)

print(
r"""
  Experiment 104 revealed terminal zero patterns in the centered parity
  coefficients.

  This experiment tests whether those zeros are explained exactly by
  a second falling-factorial law in the channel index:

      coefficient(k)
        = (K-k)_(s) * Q(k),

  where K is the final observed channel index.

  This is a genuinely different structure from the already-established

      P_k(j)
        = j_(k) (m-j)_(max(0,k-4)) R_k(j).

  A positive result would give a two-directional factorial organization:

      layer index:
          j_(k)

      endpoint index:
          (m-j)_(s)

      parity-index boundary:
          (K-k)_(s').

  The crucial tests are:

      1. exact terminal-zero multiplicity;
      2. exact polynomial divisibility;
      3. quotient degree reduction;
      4. quotient coefficient-matrix rank;
      5. exact reconstruction;
      6. a second-stage search for additional k-factors.

  This is still entirely finite and exact.  No extrapolation is assumed.
"""
)


# =============================================================================
# 11. FINAL EXACTNESS
# =============================================================================

print()
print("=" * 78)
print("11. FINAL EXACTNESS")
print("=" * 78)

print(
    f"  index_falling_divisibility = {all_divisions}"
)

print(
    f"  quotient_reconstruction = {reconstruction_ok}"
)

ok = (
    all_divisions
    and reconstruction_ok
)

print(
    f"  failures = {0 if ok else 1}"
)

print(
    f"  ALL BASIC CHECKS PASS = {ok}"
)

print()
print("EXPERIMENT 105 COMPLETE")

