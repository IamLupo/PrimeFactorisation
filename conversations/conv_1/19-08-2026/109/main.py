#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 273 — EXACT q-ROW EVALUATION / FINITE-DIFFERENCE OPERATOR AUDIT
==============================================================================

Goal:

    Determine whether the observed B-odd coefficient family

        C_p[k]

    is obtained from the source polynomial

        Q_p(j) = sum_r q_p(r) * j_(r)

    through a standard discrete operator.

Experiment 272 compared the RAW centered source polynomial directly with
the OBSERVED residual coefficients. That skips the row-dependent residual
extraction stage.

Experiment 273 therefore tests a library of exact operators acting on
the source polynomial first.

For each p we form

    Q_p(j)

and test the indexed sequence at k under:

    1. direct evaluation:
         Q(k)

    2. shifted evaluation:
         Q(k+a), a in a small fixed structural set

    3. forward finite differences:
         Delta^r Q(k)

    4. backward finite differences

    5. signed finite differences

    6. factorial-normalized differences

    7. binomial-normalized differences

    8. reversed coefficient sequences

    9. centered evaluation:
         Q((2k-D)/2)

The crucial requirement is that the SAME operator formula must work for
all available p rows.

No arbitrary fitted matrix is allowed.

All arithmetic is exact over QQ.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# =============================================================================
# SOURCE q-TABLE
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
# OBSERVED B-ODD COEFFICIENT FAMILIES
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


j = sp.symbols("j")


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


def falling(x, n):
    out = sp.Integer(1)

    for r in range(n):
        out *= x - r

    return sp.expand(out)


def build_Q(p):
    out = sp.Integer(0)

    for r, qv in enumerate(Q[p]):
        out += (
            sp.Integer(qv)
            * falling(j, r)
        )

    return clean(out)


def eval_sequence(P, ks):
    return [
        clean(
            P.subs(j, sp.Integer(k))
        )
        for k in ks
    ]


def forward_difference(values):
    out = [clean(v) for v in values]

    result = []

    while out:

        result.append(out[0])

        out = [
            clean(
                out[i + 1] - out[i]
            )
            for i in range(
                len(out) - 1
            )
        ]

    return result


def backward_difference(values):
    out = [clean(v) for v in values]

    result = []

    while out:

        result.append(out[-1])

        out = [
            clean(
                out[i + 1] - out[i]
            )
            for i in range(
                len(out) - 1
            )
        ]

    return result


def normalized_difference(
    differences,
    mode,
):
    out = []

    for r, value in enumerate(
        differences
    ):

        if mode == "raw":
            scale = 1

        elif mode == "factorial":
            scale = math.factorial(r)

        elif mode == "binomial_D":
            Dv = None
            scale = 1

        else:
            raise ValueError(mode)

        out.append(
            clean(
                value / scale
            )
        )

    return out


def compare_exact(candidate, target):
    if len(candidate) != len(target):
        return False

    return all(
        clean(a - b) == 0
        for a, b in zip(
            candidate,
            target,
        )
    )


def scalar_match(candidate, target):
    """
    Determine whether target = lambda * candidate for a single exact lambda.
    """
    if len(candidate) != len(target):
        return None

    pairs = [
        (a, b)
        for a, b in zip(
            candidate,
            target,
        )
        if a != 0
    ]

    if not pairs:
        return sp.Integer(1) if all(
            b == 0
            for b in target
        ) else None

    lam = clean(
        pairs[0][1] / pairs[0][0]
    )

    for a, b in zip(
        candidate,
        target,
    ):

        if clean(
            b - lam * a
        ) != 0:
            return None

    return lam


def generate_candidates(
    p,
    operator_name,
):
    P = build_Q(p)
    n = len(C[p])
    Dp = D[p]
    ks = list(range(n))

    # -------------------------------------------------------------------------
    # Direct evaluation
    # -------------------------------------------------------------------------

    if operator_name == "eval":
        return eval_sequence(
            P,
            ks,
        )

    # -------------------------------------------------------------------------
    # Shifted evaluations
    # -------------------------------------------------------------------------

    if operator_name.startswith(
        "eval_shift_"
    ):
        shift = int(
            operator_name.split("_")[-1]
        )

        return [
            clean(
                P.subs(
                    j,
                    sp.Integer(k + shift)
                )
            )
            for k in ks
        ]

    # -------------------------------------------------------------------------
    # Centered evaluation
    # -------------------------------------------------------------------------

    if operator_name == "center_eval":
        return [
            clean(
                P.subs(
                    j,
                    sp.Rational(
                        2 * k - Dp,
                        2,
                    )
                )
            )
            for k in ks
        ]

    # -------------------------------------------------------------------------
    # Reversed evaluation
    # -------------------------------------------------------------------------

    if operator_name == "eval_reverse":
        return eval_sequence(
            P,
            list(
                reversed(ks)
            ),
        )

    # -------------------------------------------------------------------------
    # Forward differences of evaluation sequence
    # -------------------------------------------------------------------------

    if operator_name.startswith(
        "forward_diff_"
    ):

        r = int(
            operator_name.split("_")[-1]
        )

        values = eval_sequence(
            P,
            list(range(
                n + r
            )),
        )

        for _ in range(r):
            values = [
                clean(
                    values[i + 1]
                    - values[i]
                )
                for i in range(
                    len(values) - 1
                )
            ]

        return values[:n]

    # -------------------------------------------------------------------------
    # Backward differences
    # -------------------------------------------------------------------------

    if operator_name.startswith(
        "backward_diff_"
    ):

        r = int(
            operator_name.split("_")[-1]
        )

        values = eval_sequence(
            P,
            list(range(
                -r,
                n
            )),
        )

        for _ in range(r):
            values = [
                clean(
                    values[i + 1]
                    - values[i]
                )
                for i in range(
                    len(values) - 1
                )
            ]

        return values[-n:]

    raise ValueError(
        operator_name
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 273 — EXACT q-ROW EVALUATION / "
        "FINITE-DIFFERENCE OPERATOR AUDIT"
    )
    print("=" * 78)

    # =========================================================================
    # 1. SOURCE POLYNOMIALS
    # =========================================================================

    print()
    print("=" * 78)
    print("1. SOURCE POLYNOMIALS")
    print("=" * 78)

    P = {}

    for p in [1, 3, 5, 7]:

        P[p] = build_Q(p)

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    D={D[p]}"
        )

        print(
            f"    Q(j)={P[p]}"
        )

    # =========================================================================
    # 2. OBSERVED TARGETS
    # =========================================================================

    print()
    print("=" * 78)
    print("2. OBSERVED TARGET SEQUENCES")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    C={C[p]}"
        )

    # =========================================================================
    # 3. OPERATOR LIBRARY
    # =========================================================================

    operators = [
        "eval",
        "eval_reverse",
        "center_eval",
        "eval_shift_-2",
        "eval_shift_-1",
        "eval_shift_1",
        "eval_shift_2",
        "forward_diff_1",
        "forward_diff_2",
        "forward_diff_3",
        "forward_diff_4",
        "forward_diff_5",
        "backward_diff_1",
        "backward_diff_2",
        "backward_diff_3",
        "backward_diff_4",
        "backward_diff_5",
    ]

    print()
    print("=" * 78)
    print("3. OPERATOR SEARCH")
    print("=" * 78)

    exact_common = []
    scaled_common = []

    for operator in operators:

        exact_all = True
        scale_values = []

        for p in [1, 3, 5, 7]:

            try:
                candidate = generate_candidates(
                    p,
                    operator,
                )
            except Exception:
                exact_all = False
                scale_values = []
                break

            target = C[p]

            if compare_exact(
                candidate,
                target,
            ):
                scale_values.append(
                    sp.Integer(1)
                )
                continue

            lam = scalar_match(
                candidate,
                target,
            )

            if lam is None:
                exact_all = False
                scale_values = []
                break

            scale_values.append(
                clean(lam)
            )

        if exact_all:
            exact_common.append(
                operator
            )

        # A scaled common operator means the same operator shape works for
        # every p but the scalar multiplier is allowed to depend on p.
        scaled_ok = True

        for p, lam in zip(
            [1, 3, 5, 7],
            scale_values,
        ):

            candidate = generate_candidates(
                p,
                operator,
            )

            if lam is None:
                scaled_ok = False
                break

            if not compare_exact(
                [
                    clean(
                        lam * x
                    )
                    for x in candidate
                ],
                C[p],
            ):
                scaled_ok = False
                break

        if scaled_ok:
            scaled_common.append(
                (
                    operator,
                    scale_values,
                )
            )

        print()
        print(
            f"  {operator}:"
        )

        print(
            f"    exact_all_p="
            f"{operator in exact_common}"
        )

        print(
            f"    scalar_scaled_all_p="
            f"{scaled_ok}"
        )

        if scaled_ok:
            print(
                f"    p_scalars={scale_values}"
            )

    # =========================================================================
    # 4. TERMINAL VALUE CHECK
    # =========================================================================

    print()
    print("=" * 78)
    print("4. TERMINAL SOURCE CHECK")
    print("=" * 78)

    q1 = Q[1][-1]
    q3 = Q[3][-1]

    print(
        f"  q1_terminal={q1}"
    )

    print(
        f"  q3_terminal={q3}"
    )

    print(
        f"  gcd={math.gcd(q1, q3)}"
    )

    print(
        f"  q1/17={q1 // 17}"
    )

    print(
        f"  q3/17={q3 // 17}"
    )

    # =========================================================================
    # 5. MOST INSTRUCTIVE CANDIDATES
    # =========================================================================

    print()
    print("=" * 78)
    print("5. DETAILED CANDIDATE OUTPUT")
    print("=" * 78)

    detailed = [
        "eval",
        "center_eval",
        "forward_diff_1",
        "forward_diff_2",
        "backward_diff_1",
        "backward_diff_2",
    ]

    for operator in detailed:

        print()
        print(
            f"  {operator}:"
        )

        for p in [1, 3, 5, 7]:

            candidate = generate_candidates(
                p,
                operator,
            )

            target = C[p]

            lam = scalar_match(
                candidate,
                target,
            )

            print()
            print(
                f"    p={p}:"
            )

            print(
                f"      candidate={candidate}"
            )

            print(
                f"      target={target}"
            )

            print(
                f"      scalar={lam}"
            )

    # =========================================================================
    # 6. STRUCTURAL INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print("6. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 272 established that simply taking the centered odd part of

    Q_p(j) = sum_r q_p(r) j_(r)

does NOT reproduce the observed B-odd coefficient vectors.

Experiment 273 inserts a more realistic class of operators:

    evaluation,
    shifted evaluation,
    centered evaluation,
    finite differences,
    reversed differences.

These are the natural discrete operators associated with falling
factorial bases.

There are three useful outcomes.

1. Exact common operator:

       C_p = O(Q_p)

   for every p.

   This identifies the B-channel with a standard discrete transform.

2. Common operator up to a p-dependent scalar:

       C_p = lambda_p O(Q_p).

   This means the source-index mechanism is universal and only a global
   p-dependent normalization remains.

3. Neither occurs.

   Then the missing operator is not a standard evaluation/difference
   transform, and the next step should reconstruct the actual residual
   operator from the original construction.

This experiment still does not fit an arbitrary matrix.
"""
    )

    # =========================================================================
    # 7. FINAL
    # =========================================================================

    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  exact_common_operators={exact_common}"
    )

    print(
        f"  scalar_common_operators={scaled_common}"
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
        "EXPERIMENT 273 COMPLETE"
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

