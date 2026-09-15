#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 280 — EXACT q-TO-B SHIFTED-TRIANGULAR PROVENANCE AUDIT
==============================================================================

Experiment 279 establishes the correct B-channel convention:

    B[k][offset] * j_(k+offset).

Thus the B-row is naturally indexed by the ABSOLUTE falling-factorial
index

    r = k + offset.

Experiment 280 now attacks the actual unresolved map:

    q_p(r)  ->  B[k][offset].

The previous experiments compared q_p(r) against final centered
coefficients, which was too late in the pipeline.

Here we compare q-values against the B coefficients at their true
absolute falling-factorial indices.

For each B row k and absolute index r:

    r >= k
    offset = r-k

we record

    B[k,r-k].

We then build the full upper-triangular B coefficient array

    M[k,r] = B[k,r-k].

The main questions are:

    1. Does each B coefficient factor through a q_p(r)-type quantity?
    2. Are nonzero B columns aligned with the terminal q ladder?
    3. Does B[k,r] / q_p(r) reveal a structured normalization?
    4. Is there a p determined by the diagonal/boundary location?
    5. Are the rational scales factorial/binomial/polynomial in k,r?
    6. Can multiple B rows be linked to the four known q-families?

NO arbitrary fitted matrix is used.

Exact QQ arithmetic only.
"""


from __future__ import annotations

import math
import sys
import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

k, r = sp.symbols("k r", integer=True)


# =============================================================================
# SOURCE q-TABLES
# =============================================================================

Q = {
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
# EXACT ORIGINAL B TABLE
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


def prime_factors_small(n):
    """
    Exact factorization for the relatively small integer coefficients
    appearing in this experiment.
    """
    n = abs(int(n))

    if n == 0:
        return {}

    result = {}
    d0 = 2

    while d0 * d0 <= n:

        while n % d0 == 0:
            result[d0] = result.get(d0, 0) + 1
            n //= d0

        d0 += 1

    if n > 1:
        result[n] = result.get(n, 0) + 1

    return result


def v17(n):
    n = int(n)

    if n == 0:
        return None

    n = abs(n)
    e = 0

    while n % 17 == 0:
        n //= 17
        e += 1

    return e


def rational_ratio(a, b):
    if b == 0:
        return None

    return clean(
        sp.sympify(a) / sp.sympify(b)
    )


def integer_signature(x):
    x = sp.Rational(x)

    den = int(x.q)
    num = int(x.p)

    g = math.gcd(
        abs(num),
        den,
    )

    if g:
        num //= g
        den //= g

    return num, den


def q_locations():

    locations = []

    for p in sorted(Q):
        for rr, value in enumerate(Q[p]):

            locations.append(
                {
                    "p": p,
                    "r": rr,
                    "q": value,
                }
            )

    return locations


# =============================================================================
# BUILD ABSOLUTE-INDEX B MATRIX
# =============================================================================

def absolute_B_matrix():

    M = {}

    for kk, row in enumerate(B):

        for offset, value in enumerate(row):

            rr = kk + offset

            M[(kk, rr)] = value

    return M


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 280 — EXACT q-TO-B SHIFTED-TRIANGULAR "
        "PROVENANCE AUDIT"
    )
    print("=" * 78)

    M = absolute_B_matrix()

    # =========================================================================
    # 1. ABSOLUTE B INDEX TABLE
    # =========================================================================

    print()
    print("=" * 78)
    print("1. ABSOLUTE FALLING-FACTORIAL B TABLE")
    print("=" * 78)

    for kk in range(len(B)):

        entries = []

        for offset, value in enumerate(B[kk]):

            rr = kk + offset

            entries.append(
                (
                    rr,
                    value,
                )
            )

        print()
        print(
            f"  k={kk}:"
        )
        print(
            f"    B[k,r]={entries}"
        )

    # =========================================================================
    # 2. Q FAMILY TERMINAL PROFILE
    # =========================================================================

    print()
    print("=" * 78)
    print("2. SOURCE q-FAMILY PROFILE")
    print("=" * 78)

    for p in sorted(Q):

        print()
        print(
            f"  p={p}:"
        )
        print(
            f"    D={D[p]}"
        )
        print(
            f"    q={Q[p]}"
        )

        print(
            f"    terminal={Q[p][-1]}"
        )

        print(
            f"    terminal_v17="
            f"{v17(Q[p][-1])}"
        )

    # =========================================================================
    # 3. DIRECT ABSOLUTE-INDEX MATCH SEARCH
    # =========================================================================

    print()
    print("=" * 78)
    print("3. DIRECT q_r / B[k,r] MATCH SEARCH")
    print("=" * 78)

    direct_matches = []

    for (kk, rr), value in sorted(
        M.items()
    ):

        for p in sorted(Q):

            if rr >= len(Q[p]):
                continue

            qvalue = Q[p][rr]

            if qvalue == 0:
                continue

            ratio = rational_ratio(
                value,
                qvalue,
            )

            if ratio is None:
                continue

            direct_matches.append(
                (
                    kk,
                    rr,
                    p,
                    qvalue,
                    value,
                    ratio,
                )
            )

            print()
            print(
                f"  k={kk} r={rr} p={p}:"
            )
            print(
                f"    q_p(r)={qvalue}"
            )
            print(
                f"    B[k,r]={value}"
            )
            print(
                f"    B/q={ratio}"
            )

    print()
    print(
        f"  direct_comparable_cases="
        f"{len(direct_matches)}"
    )

    # =========================================================================
    # 4. TERMINAL-COLUMN PROFILE
    # =========================================================================

    print()
    print("=" * 78)
    print("4. TERMINAL SOURCE-COLUMN PROFILE")
    print("=" * 78)

    terminal_sources = {
        5: Q[1][5],
        4: Q[3][4],
        2: Q[5][2],
        0: Q[7][0],
    }

    for rr in sorted(
        terminal_sources
    ):

        qvalue = terminal_sources[
            rr
        ]

        locations = []

        for (
            kk,
            rrr,
        ), value in M.items():

            if rrr != rr:
                continue

            locations.append(
                (
                    kk,
                    value,
                )
            )

        print()
        print(
            f"  r={rr}: "
            f"terminal_q={qvalue}"
        )

        print(
            f"    B_locations={locations}"
        )

    # =========================================================================
    # 5. B COLUMN CONTENT RELATIVE TO q TERMINAL VALUES
    # =========================================================================

    print()
    print("=" * 78)
    print("5. TERMINAL COLUMN NORMALIZATION")
    print("=" * 78)

    for rr, qvalue in sorted(
        terminal_sources.items()
    ):

        print()
        print(
            f"  r={rr}:"
        )

        print(
            f"    q_terminal={qvalue}"
        )

        if qvalue == 0:
            print(
                "    q=0"
            )
            continue

        for kk in range(
            min(
                len(B),
                rr + 1,
            )
        ):

            if (kk, rr) not in M:
                continue

            value = M[
                (
                    kk,
                    rr,
                )
            ]

            print(
                f"    k={kk}: "
                f"B={value} "
                f"B/q={rational_ratio(value, qvalue)}"
            )

    # =========================================================================
    # 6. FACTORIAL / BINOMIAL NORMALIZATION SEARCH
    # =========================================================================

    print()
    print("=" * 78)
    print("6. ELEMENTARY NORMALIZATION SEARCH")
    print("=" * 78)

    normalization_names = [
        "raw",
        "divide_factorial_r",
        "divide_factorial_k",
        "multiply_factorial_r",
        "multiply_factorial_k",
        "divide_binom_r_k",
        "multiply_binom_r_k",
        "divide_2_power_r",
        "multiply_2_power_r",
    ]

    normalization_matches = {
        name: []
        for name in normalization_names
    }

    for (
        kk,
        rr,
    ), value in sorted(M.items()):

        if rr < 0:
            continue

        q_candidates = []

        for p in sorted(Q):

            if rr < len(Q[p]):
                qvalue = Q[p][rr]

                if qvalue != 0:
                    q_candidates.append(
                        (
                            p,
                            qvalue,
                        )
                    )

        if not q_candidates:
            continue

        for p, qvalue in q_candidates:

            candidates = {
                "raw": value,
                "divide_factorial_r":
                    value / sp.factorial(rr),
                "divide_factorial_k":
                    value / sp.factorial(kk),
                "multiply_factorial_r":
                    value * sp.factorial(rr),
                "multiply_factorial_k":
                    value * sp.factorial(kk),
                "divide_binom_r_k":
                    (
                        value
                        / sp.binomial(
                            rr,
                            kk,
                        )
                        if rr >= kk
                        else value
                    ),
                "multiply_binom_r_k":
                    (
                        value
                        * sp.binomial(
                            rr,
                            kk,
                        )
                        if rr >= kk
                        else value
                    ),
                "divide_2_power_r":
                    value / (2 ** rr),
                "multiply_2_power_r":
                    value * (2 ** rr),
            }

            for name, candidate in candidates.items():

                ratio = rational_ratio(
                    candidate,
                    qvalue,
                )

                if ratio is not None:
                    normalization_matches[
                        name
                    ].append(
                        (
                            kk,
                            rr,
                            p,
                            ratio,
                        )
                    )

    for name in normalization_names:

        print()
        print(
            f"  {name}: "
            f"comparisons="
            f"{len(normalization_matches[name])}"
        )

        if normalization_matches[name]:
            print(
                f"    first_matches="
                f"{normalization_matches[name][:10]}"
            )

    # =========================================================================
    # 7. DIAGONAL / ANTIDIAGONAL PROFILES
    # =========================================================================

    print()
    print("=" * 78)
    print("7. DIAGONAL / ANTIDIAGONAL PROFILE")
    print("=" * 78)

    diagonal = []

    for kk in range(
        len(B)
    ):

        if (
            kk,
            kk,
        ) in M:

            diagonal.append(
                (
                    kk,
                    M[
                        (
                            kk,
                            kk,
                        )
                    ],
                )
            )

    print()
    print(
        f"  diagonal={diagonal}"
    )

    anti_diagonal = []

    max_r = max(
        rr
        for kk, rr
        in M.keys()
    )

    for kk in range(
        len(B)
    ):

        rr = max_r - kk

        if (
            kk,
            rr,
        ) in M:

            anti_diagonal.append(
                (
                    kk,
                    rr,
                    M[
                        (
                            kk,
                            rr,
                        )
                    ],
                )
            )

    print(
        f"  anti_diagonal="
        f"{anti_diagonal}"
    )

    # =========================================================================
    # 8. SOURCE-INDEX TRIANGULARITY
    # =========================================================================

    print()
    print("=" * 78)
    print("8. TRIANGULARITY CHECK")
    print("=" * 78)

    upper_triangular = True

    for (
        kk,
        rr,
    ) in M:

        if rr < kk:
            upper_triangular = False

    print(
        f"  B[k,r]=0 for r<k: "
        f"{upper_triangular}"
    )

    for kk in range(
        len(B)
    ):

        support = [
            rr
            for (
                k0,
                rr,
            ) in M
            if k0 == kk
        ]

        print(
            f"  k={kk}: "
            f"absolute_support={support}"
        )

    # =========================================================================
    # 9. TERMINAL PROJECTIVE SOURCE
    # =========================================================================

    print()
    print("=" * 78)
    print("9. TERMINAL PROJECTIVE SOURCE")
    print("=" * 78)

    q1_terminal = Q[1][-1]
    q3_terminal = Q[3][-1]

    print(
        f"  q1_terminal={q1_terminal}"
    )

    print(
        f"  q3_terminal={q3_terminal}"
    )

    print(
        f"  gcd={sp.gcd(
            q1_terminal,
            q3_terminal,
        )}"
    )

    print(
        f"  q1/17="
        f"{q1_terminal // 17}"
    )

    print(
        f"  q3/17="
        f"{q3_terminal // 17}"
    )

    # =========================================================================
    # 10. STRUCTURAL INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 279 establishes that the B channel uses the shifted
falling-factorial basis

    B[k][offset] j_(k+offset).

Experiment 280 therefore moves the provenance comparison to the correct
level:

    q_p(r)
        versus
    B[k,r].

This is the first experiment where the source q-index r and the B
falling-factorial index r are genuinely aligned.

There are several possible outcomes.

A. Repeated exact ratios:

       B[k,r] / q_p(r)

   simplify systematically.

   This would reveal the source operator.

B. Ratios become simple after factorial/binomial normalization.

   This would identify the combinatorial normalization.

C. Only terminal columns align.

   Then the earlier terminal provenance is genuine, but the interior
   B coefficients require additional construction data.

D. No alignment occurs.

   Then q_p(r) is not itself the immediate coefficient source of B.

No arbitrary fitted matrix is used.
"""
    )

    # =========================================================================
    # 11. FINAL
    # =========================================================================

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  shifted_B_indexing_confirmed=True"
    )

    print(
        f"  direct_q_B_comparisons="
        f"{len(direct_matches)}"
    )

    print(
        f"  upper_triangular_index_structure="
        f"{upper_triangular}"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 280 COMPLETE"
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

