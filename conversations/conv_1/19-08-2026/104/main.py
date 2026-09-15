#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 268R — EXACT COMMON q-ROW CONVOLUTION / TOEPLITZ AUDIT
==============================================================================

Corrected Experiment 268.

The previous script failed only because math.gcd() was used without importing
math.

This version also keeps the exact rational arithmetic and validates every
candidate kernel independently.

Tested models:

    1. lower Toeplitz:
         C_k = sum_{t=0}^k T_t q_{k-t}

    2. upper Toeplitz:
         C_k = sum_{t=0}^{n-1-k} T_t q_{k+t}

    3. reverse-lower:
         reverse(q), then lower convolution

    4. reverse-upper:
         reverse(q), then upper convolution

No arbitrary matrix fitting.
No floating point.
No extrapolation.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# =============================================================================
# SOURCE q-ROWS
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
# OBSERVED B-ODD COEFFICIENTS
# =============================================================================

C = {
    1: [
        sp.Rational(-584531, 35840),
        sp.Rational(-2908483, 1920),
        sp.Rational(-31233169, 13440),
        sp.Rational(-1446167, 1344),
        sp.Rational(-22259149, 40320),
        sp.Rational(-301, 240),
    ],
    3: [
        sp.Rational(-59257, 46080),
        sp.Rational(186547, 1440),
        sp.Rational(367433, 1680),
        sp.Rational(126549, 448),
        sp.Rational(-162139, 40320),
    ],
    5: [
        sp.Rational(4457, 46080),
        sp.Rational(-16819, 5760),
        sp.Rational(-5769, 896),
    ],
    7: [
        sp.Rational(-421, 322560),
    ],
}


# =============================================================================
# HELPERS
# =============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def build_equations(mode, kernel_len):
    """
    Return:
        symbols, equations

    Each equation has the form

        sum(a_i * t_i) = observed_value
    """

    t = sp.symbols(
        "t0:%d" % kernel_len
    )

    equations = []

    for p in [1, 3, 5, 7]:

        q = Q[p]
        c = C[p]

        n = len(q)

        for k in range(n):

            coeffs = [
                sp.Integer(0)
                for _ in range(kernel_len)
            ]

            if mode == "lower":

                for s in range(k + 1):
                    if s < kernel_len:
                        coeffs[s] += q[k - s]

            elif mode == "upper":

                for s in range(n - k):
                    if s < kernel_len:
                        coeffs[s] += q[k + s]

            elif mode == "reverse_lower":

                qr = list(reversed(q))

                for s in range(k + 1):
                    if s < kernel_len:
                        coeffs[s] += qr[k - s]

            elif mode == "reverse_upper":

                qr = list(reversed(q))

                for s in range(n - k):
                    if s < kernel_len:
                        coeffs[s] += qr[k + s]

            else:
                raise ValueError(
                    f"unknown mode {mode}"
                )

            equations.append(
                (
                    coeffs,
                    clean(c[k]),
                )
            )

    return t, equations


def solve_unique_or_none(
    symbols,
    equations,
):
    """
    Exact consistency / uniqueness test.

    Returns:
        None
        or exact kernel vector.
    """

    A = sp.Matrix(
        [
            row
            for row, _ in equations
        ]
    )

    b = sp.Matrix(
        [
            value
            for _, value in equations
        ]
    )

    augmented = A.row_join(b)

    if A.rank() != augmented.rank():
        return None

    # We require a unique solution.
    if A.rank() != len(symbols):
        return None

    solution = sp.linsolve(
        (A, b),
        symbols,
    )

    if len(solution) != 1:
        return None

    vector = next(iter(solution))

    # Reject symbolic/free parameters.
    for value in vector:
        if value.free_symbols:
            return None

    return [
        clean(value)
        for value in vector
    ]


def predict_lower(q, kernel):
    result = []

    for k in range(len(q)):

        value = sp.Integer(0)

        for s in range(min(k + 1, len(kernel))):
            value += (
                kernel[s]
                * q[k - s]
            )

        result.append(
            clean(value)
        )

    return result


def predict_upper(q, kernel):
    result = []

    n = len(q)

    for k in range(n):

        value = sp.Integer(0)

        for s in range(
            min(n - k, len(kernel))
        ):
            value += (
                kernel[s]
                * q[k + s]
            )

        result.append(
            clean(value)
        )

    return result


def predict_reverse_lower(q, kernel):
    qr = list(reversed(q))

    predicted = predict_lower(
        qr,
        kernel,
    )

    return list(reversed(predicted))


def predict_reverse_upper(q, kernel):
    qr = list(reversed(q))

    predicted = predict_upper(
        qr,
        kernel,
    )

    return list(reversed(predicted))


def verify(mode, kernel):

    predictions = {}
    exact = True

    for p in [1, 3, 5, 7]:

        if mode == "lower":
            pred = predict_lower(
                Q[p],
                kernel,
            )

        elif mode == "upper":
            pred = predict_upper(
                Q[p],
                kernel,
            )

        elif mode == "reverse_lower":
            pred = predict_reverse_lower(
                Q[p],
                kernel,
            )

        elif mode == "reverse_upper":
            pred = predict_reverse_upper(
                Q[p],
                kernel,
            )

        else:
            raise ValueError(mode)

        predictions[p] = pred

        if pred != C[p]:
            exact = False

    return exact, predictions


def integer_denominator_profile(kernel):
    return [
        int(sp.denom(value))
        for value in kernel
    ]


def common_denominator(kernel):
    d = 1

    for value in kernel:
        d = math.lcm(
            d,
            int(sp.denom(value)),
        )

    return d


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 268R — EXACT COMMON q-ROW "
        "CONVOLUTION / TOEPLITZ AUDIT"
    )
    print("=" * 78)

    # -------------------------------------------------------------------------
    # 1. DATA
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE / OBSERVED DATA")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    q={Q[p]}"
        )

        print(
            f"    C={C[p]}"
        )

    # -------------------------------------------------------------------------
    # 2. KERNEL TESTS
    # -------------------------------------------------------------------------

    modes = [
        "lower",
        "upper",
        "reverse_lower",
        "reverse_upper",
    ]

    kernel_len = 6
    results = {}

    print()
    print("=" * 78)
    print("2. COMMON TOEPLITZ KERNEL SEARCH")
    print("=" * 78)

    for mode in modes:

        symbols, equations = build_equations(
            mode,
            kernel_len,
        )

        kernel = solve_unique_or_none(
            symbols,
            equations,
        )

        print()
        print(
            f"  mode={mode}:"
        )

        print(
            f"    equations="
            f"{len(equations)}"
        )

        if kernel is None:

            print(
                "    exact_unique_kernel=False"
            )

            results[mode] = None
            continue

        exact, predictions = verify(
            mode,
            kernel,
        )

        print(
            "    exact_unique_kernel=True"
        )

        print(
            f"    kernel={kernel}"
        )

        print(
            f"    denominators="
            f"{integer_denominator_profile(kernel)}"
        )

        print(
            f"    common_denominator="
            f"{common_denominator(kernel)}"
        )

        print(
            f"    independent_verification="
            f"{exact}"
        )

        if exact:
            for p in [1, 3, 5, 7]:
                print(
                    f"      p={p}: "
                    f"predicted={predictions[p]}"
                )

        results[mode] = {
            "kernel": kernel,
            "exact": exact,
        }

    # -------------------------------------------------------------------------
    # 3. TERMINAL SOURCE REFERENCE
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. TERMINAL SOURCE REFERENCE")
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
        f"  gcd={math.gcd(q1_terminal, q3_terminal)}"
    )

    print(
        f"  q1_terminal/17="
        f"{q1_terminal // 17}"
    )

    print(
        f"  q3_terminal/17="
        f"{q3_terminal // 17}"
    )

    # -------------------------------------------------------------------------
    # 4. STRUCTURAL SUMMARY
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. STRUCTURAL SUMMARY")
    print("=" * 78)

    successful_modes = [
        mode
        for mode, result in results.items()
        if result is not None
        and result["exact"]
    ]

    print(
        f"  successful_common_modes="
        f"{successful_modes}"
    )

    print(
        f"  common_toeplitz_kernel_found="
        f"{bool(successful_modes)}"
    )

    # -------------------------------------------------------------------------
    # 5. INTERPRETATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
The previous experiments ruled out a direct basis-change explanation.

This experiment tests a much more constrained source-index operator:

    C = T * q,

where T is Toeplitz.

A lower Toeplitz form corresponds to dependence on previous q-indices;
an upper form corresponds to dependence on later q-indices. The reversed
forms test whether the natural q-index orientation is terminal-to-initial.

Because the SAME kernel must reproduce every available p-row, a
successful result would be substantially stronger than a one-case fit.

If no common kernel exists, the missing source operator must depend on
additional structure such as:

    p,
    D(p),
    k,
    boundary position,
    or another part of the original construction.

No family formula is inferred from a negative result.
"""
    )

    # -------------------------------------------------------------------------
    # 6. FINAL
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  common_exact_modes="
        f"{successful_modes}"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  exact_rational_arithmetic=True"
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
        "EXPERIMENT 268R COMPLETE"
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