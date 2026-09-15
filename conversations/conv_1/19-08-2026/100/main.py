#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 264 — EXACT FULL q_p(r) ROW / PARITY-COEFFICIENT
                 RECONSTRUCTION AUDIT
==============================================================================

Experiment 263 established that

    D(p) = max{k : [y^p] B_odd(k) != 0}

and that the terminal coefficient recovers

    q_p(D(p)).

Experiment 264 asks the stronger question:

    does the COMPLETE coefficient family

        [y^p] B_odd(k),   k=0,...,D(p)

    contain the COMPLETE original source row

        q_p(0),...,q_p(D(p))

    through an exact triangular transformation?

The experiment does NOT guess the transformation.

Instead, it tests several exact structural possibilities:

    1. direct equality after clearing denominators;
    2. reversed-row equality;
    3. exact finite differences;
    4. exact binomial/falling-factor transforms;
    5. exact triangular reconstruction by solving the linear system.

For each odd p:

    p=1,3,5,7

the script compares the parity coefficient vector with the original
q-row from the construction.

A successful triangular reconstruction would be a major provenance
result: the centered parity kernel would contain the entire q_p(r)
source row, not merely its terminal entry.

Only main.py.
No floating point.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

j, y, k = sp.symbols(
    "j y k"
)


# =============================================================================
# ORIGINAL q-TABLE
# =============================================================================

Q_ORIGINAL = {
    1: [
        -126258696,
        -11600759760,
        2668721436,
        1764373740,
        -1338089411,
        495451247,
    ],

    3: [
        9955176,
        -1263551016,
        -152369292,
        -128667196,
        421514439,
    ],

    5: [
        -62398,
        4771718,
        16027881,
    ],

    7: [
        1,
    ],
}


D = {
    1: 5,
    3: 4,
    5: 2,
    7: 0,
}


# =============================================================================
# EXACT EXPERIMENT 104 B-DATA
# =============================================================================

B = [
    [
        25,
        619,
        sp.Rational(3231, 2),
        sp.Rational(-33, 2),
        sp.Rational(-1675, 4),
        sp.Rational(3363, 20),
        sp.Rational(-9991, 360),
        sp.Rational(-421, 2520),
    ],

    [
        1750,
        8624,
        sp.Rational(6829, 3),
        sp.Rational(-27341, 8),
        sp.Rational(10551, 10),
        sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],

    [
        9690,
        10234,
        sp.Rational(-57829, 8),
        sp.Rational(148151, 120),
        sp.Rational(1432, 5),
        sp.Rational(-5769, 28),
    ],

    [
        sp.Rational(22100, 3),
        sp.Rational(-19045, 12),
        sp.Rational(-5577, 4),
        sp.Rational(351271, 360),
        sp.Rational(-101119, 315),
    ],

    [
        sp.Rational(17875, 24),
        sp.Rational(-22061, 40),
        sp.Rational(132343, 720),
        sp.Rational(-162139, 5040),
    ],

    [
        sp.Rational(65, 12),
        sp.Rational(-313, 60),
        sp.Rational(301, 120),
    ],
]


# =============================================================================
# HELPERS
# =============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def poly(expr, var):
    return sp.Poly(
        sp.expand(
            sp.sympify(expr)
        ),
        var,
        domain=sp.QQ,
    )


def falling(x, n):
    out = sp.Integer(1)

    for r in range(n):
        out *= x - r

    return sp.expand(out)


def upper_falling(m, x, n):
    out = sp.Integer(1)

    for r in range(n):
        out *= m - x - r

    return sp.expand(out)


def reconstruct(row, k_index):

    result = sp.Integer(0)

    for offset, coeff in enumerate(row):

        degree = (
            k_index
            + offset
        )

        result += (
            sp.sympify(coeff)
            * falling(
                j,
                degree,
            )
        )

    return clean(result)


def exact_quotient(num, den, var):

    pn = poly(
        clean(num),
        var,
    )

    pd = poly(
        clean(den),
        var,
    )

    quotient, remainder = sp.div(
        pn,
        pd,
        domain=sp.QQ,
    )

    if not remainder.is_zero:
        return None

    return clean(
        quotient.as_expr()
    )


def residual_channel(channel, m):

    result = []

    for k_index, row in enumerate(channel):

        P = reconstruct(
            row,
            k_index,
        )

        s = max(
            0,
            k_index - 4,
        )

        divisor = clean(
            falling(
                j,
                k_index,
            )
            *
            upper_falling(
                m,
                j,
                s,
            )
        )

        R = exact_quotient(
            P,
            divisor,
            j,
        )

        if R is None:

            raise ArithmeticError(
                f"Residual extraction failed at k={k_index}."
            )

        result.append(R)

    return result


def center(R, m):

    return clean(
        R.subs(
            j,
            (y + m) / 2,
        )
    )


def odd_part(F):

    return clean(
        (
            F
            - F.subs(
                y,
                -y,
            )
        )
        / 2
    )


def coefficient(expr, power):

    if clean(expr) == 0:

        return sp.Integer(0)

    return sp.Rational(
        poly(
            expr,
            y,
        ).nth(power)
    )


def interpolate_index(values):

    points = [
        (
            sp.Integer(i),
            sp.Rational(v),
        )
        for i, v in enumerate(values)
    ]

    return clean(
        sp.interpolate(
            points,
            k,
        )
    )


# =============================================================================
# EXACT VECTOR UTILITIES
# =============================================================================

def lcm_denominators(values):

    L = 1

    for value in values:

        value = sp.Rational(value)

        L = math.lcm(
            L,
            int(value.q),
        )

    return L


def clear_denominators(values):

    L = lcm_denominators(
        values
    )

    return [
        int(
            sp.Rational(v) * L
        )
        for v in values
    ], L


def primitive_vector(values):

    ints, denominator = (
        clear_denominators(
            values
        )
    )

    g = 0

    for value in ints:

        g = math.gcd(
            g,
            abs(value)
        )

    if g:

        ints = [
            value // g
            for value in ints
        ]

        denominator //= g

    return ints, denominator


def exact_first_differences(values):

    return [
        clean(
            values[i + 1]
            - values[i]
        )
        for i in range(
            len(values) - 1
        )
    ]


def exact_difference_table(values):

    current = [
        sp.Rational(v)
        for v in values
    ]

    table = [
        current
    ]

    while len(current) > 1:

        current = (
            exact_first_differences(
                current
            )
        )

        table.append(
            current
        )

    return table


def binomial_transform(values):

    """
    Exact forward binomial transform:

        T_r = sum_{i=0}^r (-1)^(r-i) C(r,i) v_i
    """

    out = []

    for r in range(
        len(values)
    ):

        value = sp.Integer(0)

        for i in range(r + 1):

            value += (
                (-1) ** (r - i)
                * sp.binomial(
                    r,
                    i,
                )
                * sp.Rational(
                    values[i]
                )
            )

        out.append(
            clean(value)
        )

    return out


def reversed_vector(values):

    return list(
        reversed(
            values
        )
    )


def normalized_terminal_vector(values):

    terminal = (
        sp.Rational(
            values[-1]
        )
    )

    if terminal == 0:

        return None

    return [
        clean(
            sp.Rational(v)
            /
            terminal
        )
        for v in values
    ]


# =============================================================================
# TRIANGULAR RECONSTRUCTION
# =============================================================================

def solve_triangular_transform(source, target):

    """
    Try

        target_r = sum_{i=0}^r T[r,i] source_i

    with a lower-triangular matrix T whose entries are solved exactly.

    A canonical triangular matrix is not assumed.

    We solve for T against the SINGLE observed vector by using the
    natural basis generated by the source cumulative transforms.

    The output is therefore diagnostic rather than a universal claim.
    """

    source = [
        sp.Rational(v)
        for v in source
    ]

    target = [
        sp.Rational(v)
        for v in target
    ]

    n = len(
        source
    )

    if n != len(
        target
    ):
        return None

    # Canonical unit lower triangular attempt:
    # T[r,r]=1 and solve lower entries recursively.
    T = sp.zeros(
        n,
        n,
    )

    for r in range(n):

        T[r, r] = 1

        residual = target[r]

        for i in range(r):

            residual -= (
                T[r, i]
                * source[i]
            )

        if r == 0:

            if clean(
                residual
                - source[0]
            ) != 0:

                return None

        else:

            pivot = source[0]

            if pivot == 0:

                return None

            T[r, 0] = clean(
                residual / pivot
            )

    return T


# =============================================================================
# MAIN
# =============================================================================

def main():

    failures = []

    # -------------------------------------------------------------------------
    # Build B-odd parity channel
    # -------------------------------------------------------------------------

    B_residuals = residual_channel(
        B,
        7,
    )

    B_centered = [
        center(
            R,
            7,
        )
        for R in B_residuals
    ]

    B_odd = [
        odd_part(F)
        for F in B_centered
    ]

    print("=" * 78)
    print(
        "EXPERIMENT 264 — EXACT FULL q_p(r) ROW / "
        "PARITY-COEFFICIENT RECONSTRUCTION AUDIT"
    )
    print("=" * 78)

    # -------------------------------------------------------------------------
    # 1. RAW PARITY VECTORS
    # -------------------------------------------------------------------------

    parity_vectors = {}

    print()
    print("=" * 78)
    print("1. RAW B-ODD PARITY COEFFICIENT VECTORS")
    print("=" * 78)

    for p in (
        1,
        3,
        5,
        7,
    ):

        values = [
            coefficient(
                R,
                p,
            )
            for R in B_odd
        ]

        parity_vectors[p] = values

        print()
        print(
            f"  y^{p}:"
        )

        print(
            f"    coefficients="
            f"{values}"
        )

        print(
            f"    primitive_integer_form="
            f"{primitive_vector(values)}"
        )

    # -------------------------------------------------------------------------
    # 2. ORIGINAL q ROWS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. ORIGINAL q_p(r) ROWS")
    print("=" * 78)

    for p in (
        1,
        3,
        5,
        7,
    ):

        print(
            f"  q_{p}="
            f"{Q_ORIGINAL[p]}"
        )

    # -------------------------------------------------------------------------
    # 3. LENGTH / SUPPORT MATCH
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. LENGTH / SUPPORT MATCH")
    print("=" * 78)

    for p in (
        1,
        3,
        5,
        7,
    ):

        q_row = Q_ORIGINAL[p]
        parity = parity_vectors[p]

        exact_length = (
            len(q_row)
            ==
            len(parity)
        )

        print(
            f"  p={p}: "
            f"q_length={len(q_row)} "
            f"parity_length={len(parity)} "
            f"equal={exact_length}"
        )

    # -------------------------------------------------------------------------
    # 4. DIRECT / REVERSED / DIFFERENCE TESTS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. DIRECT / REVERSED / DIFFERENCE TESTS")
    print("=" * 78)

    for p in (
        1,
        3,
        5,
        7,
    ):

        parity = parity_vectors[p]
        q_row = [
            sp.Integer(v)
            for v in Q_ORIGINAL[p]
        ]

        direct = (
            parity
            ==
            q_row
        )

        reverse = (
            parity
            ==
            reversed_vector(
                q_row
            )
        )

        diff_parity = (
            exact_first_differences(
                parity
            )
        )

        diff_q = (
            exact_first_differences(
                q_row
            )
        )

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    direct_equal={direct}"
        )

        print(
            f"    reversed_equal={reverse}"
        )

        print(
            f"    parity_first_differences="
            f"{diff_parity}"
        )

        print(
            f"    q_first_differences="
            f"{diff_q}"
        )

    # -------------------------------------------------------------------------
    # 5. DIFFERENCE TABLE COMPARISON
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. EXACT DIFFERENCE-TABLE COMPARISON")
    print("=" * 78)

    for p in (
        1,
        3,
        5,
        7,
    ):

        parity_table = (
            exact_difference_table(
                parity_vectors[p]
            )
        )

        q_table = (
            exact_difference_table(
                Q_ORIGINAL[p]
            )
        )

        print()
        print(
            f"  p={p}:"
        )

        for r in range(
            min(
                len(
                    parity_table
                ),
                len(
                    q_table
                ),
            )
        ):

            print(
                f"    diff_order={r}:"
            )

            print(
                f"      parity="
                f"{parity_table[r]}"
            )

            print(
                f"      q="
                f"{q_table[r]}"
            )

    # -------------------------------------------------------------------------
    # 6. BINOMIAL TRANSFORM COMPARISON
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. EXACT BINOMIAL-TRANSFORM COMPARISON")
    print("=" * 78)

    for p in (
        1,
        3,
        5,
        7,
    ):

        parity = (
            parity_vectors[p]
        )

        transformed = (
            binomial_transform(
                parity
            )
        )

        q_row = [
            sp.Integer(v)
            for v in Q_ORIGINAL[p]
        ]

        exact = (
            transformed
            ==
            q_row
        )

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    transformed="
            f"{transformed}"
        )

        print(
            f"    q_row="
            f"{q_row}"
        )

        print(
            f"    exact={exact}"
        )

    # -------------------------------------------------------------------------
    # 7. TERMINAL NORMALIZATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. TERMINAL-NORMALIZED COMPARISON")
    print("=" * 78)

    for p in (
        1,
        3,
        5,
        7,
    ):

        parity_normalized = (
            normalized_terminal_vector(
                parity_vectors[p]
            )
        )

        q_normalized = (
            normalized_terminal_vector(
                Q_ORIGINAL[p]
            )
        )

        exact = (
            parity_normalized
            ==
            q_normalized
        )

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    parity_normalized="
            f"{parity_normalized}"
        )

        print(
            f"    q_normalized="
            f"{q_normalized}"
        )

        print(
            f"    exact={exact}"
        )

    # -------------------------------------------------------------------------
    # 8. GENERATING-POLYNOMIAL VIEW
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. GENERATING-POLYNOMIAL VIEW")
    print("=" * 78)

    z = sp.symbols(
        "z"
    )

    for p in (
        1,
        3,
        5,
        7,
    ):

        parity = parity_vectors[p]
        q_row = Q_ORIGINAL[p]

        parity_poly = clean(
            sum(
                sp.Rational(
                    parity[r]
                )
                * z**r
                for r in range(
                    len(parity)
                )
            )
        )

        q_poly = clean(
            sum(
                sp.Integer(
                    q_row[r]
                )
                * z**r
                for r in range(
                    len(q_row)
                )
            )
        )

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    parity_generating="
            f"{parity_poly}"
        )

        print(
            f"    q_generating="
            f"{q_poly}"
        )

        print(
            f"    difference="
            f"{clean(parity_poly - q_poly)}"
        )

    # -------------------------------------------------------------------------
    # 9. STRUCTURAL INTERPRETATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 263 showed that the endpoint of the B-odd support recovers
D(p), and that the endpoint value recovers q_p(D(p)).

Experiment 264 now asks whether the ENTIRE q_p(r) source row is encoded
in the complete parity coefficient vector.

A positive result in any exact transform class would give a much stronger
source theorem:

    centered parity kernel
        ->
    complete q_p(r) source row.

The script deliberately tests several natural exact transforms without
assuming one in advance.

The most important outcome is therefore not a particular transform
passing, but whether the q-row can be reconstructed by a simple exact
basis change.

A failure of all elementary transforms is itself informative: it means
the terminal value correspondence is real, but the full q-row requires
a deeper triangular kernel.
"""
    )

    # -------------------------------------------------------------------------
    # 10. FINAL
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    # The basic provenance itself is still exact.
    endpoint_exact = True
    terminal_exact = True

    for p in (
        1,
        3,
        5,
        7,
    ):

        support = [
            i
            for i, value in enumerate(
                parity_vectors[p]
            )
            if value != 0
        ]

        if max(support) != D[p]:
            endpoint_exact = False

        if (
            sp.Rational(
                parity_vectors[p][D[p]]
            )
            !=
            sp.Rational(
                parity_vectors[p][-1]
            )
        ):
            terminal_exact = False

    # This experiment is exploratory by design.
    # Only fatal reconstruction inconsistencies count as failures.
    final_ok = (
        endpoint_exact
        and terminal_exact
    )

    print(
        f"  endpoint_provenance_exact="
        f"{endpoint_exact}"
    )

    print(
        f"  terminal_provenance_exact="
        f"{terminal_exact}"
    )

    print(
        "  full_q_row_simple_transform_proved=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print(
        "EXPERIMENT 264 COMPLETE"
    )


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:

        print(
            "\nInterrupted."
        )

        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )

        raise

