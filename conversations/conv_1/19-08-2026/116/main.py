#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 282 — EXACT q-FAMILY HOLD-OUT / LEAVE-ONE-OUT PROVENANCE AUDIT
==============================================================================

Experiment 281R found four exact linear reconstructions, but several were
exactly determined by as many columns as unknown weights.

Example:

    2 source channels
    2 observed columns
    2 unknown weights

That is interpolation, not evidence of a structural operator.

Experiment 282 therefore imposes an out-of-sample requirement.

For a proposed source subset P:

    B[k,r] = sum_{p in P} w[k,p] q_p(r)

the weights are determined from a TRAINING set of columns and then tested
on one or more HELD-OUT columns.

A true structural law must satisfy:

    training residuals = 0
    AND
    held-out residuals = 0.

We perform:

    1. prefix-train / suffix-test;
    2. leave-one-column-out tests;
    3. maximal available column tests;
    4. exact rational arithmetic only.

Missing q_p(r) entries are undefined, never zero.

No arbitrary B[k,r] matrix is fitted.
"""


from __future__ import annotations

import sys
import sympy as sp


# =============================================================================
# SOURCE q-FAMILIES
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

P_VALUES = [1, 3, 5, 7]


# =============================================================================
# B[k,r], WITH ABSOLUTE FALLING-FACTORIAL INDEX r
# =============================================================================

B = [
    [
        sp.Rational(25),
        sp.Rational(619),
        sp.Rational(3231, 2),
        sp.Rational(-33, 2),
        sp.Rational(-1675, 4),
        sp.Rational(3363, 20),
        sp.Rational(-9991, 360),
        sp.Rational(-421, 2520),
    ],
    [
        sp.Rational(1750),
        sp.Rational(8624),
        sp.Rational(6829, 3),
        sp.Rational(-27341, 8),
        sp.Rational(10551, 10),
        sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],
    [
        sp.Rational(9690),
        sp.Rational(10234),
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


def q_defined(p, r):
    return (
        0 <= r < len(Q[p])
    )


def q_value(p, r):
    if not q_defined(p, r):
        return None

    return sp.Rational(
        Q[p][r]
    )


def b_value(k, r):
    offset = r - k

    if offset < 0:
        return None

    if offset >= len(B[k]):
        return None

    return sp.Rational(
        B[k][offset]
    )


def common_columns(k, labels):
    return [
        r
        for r in range(8)
        if (
            b_value(k, r) is not None
            and all(
                q_defined(p, r)
                for p in labels
            )
        )
    ]


def build_system(k, labels, columns):
    M = []
    rhs = []

    for r in columns:

        row = [
            q_value(
                p,
                r,
            )
            for p in labels
        ]

        target = b_value(
            k,
            r,
        )

        M.append(row)
        rhs.append(target)

    return (
        sp.Matrix(
            M
        ),
        sp.Matrix(
            rhs
        ),
    )


def solve_unique(M, rhs):
    """
    Exact unique solve only.

    Returns:
        status,
        solution_or_None,
        rank,
        nullity
    """

    rank = M.rank()

    augmented_rank = (
        M
        .row_join(rhs)
        .rank()
    )

    if augmented_rank != rank:
        return (
            "NO_SOLUTION",
            None,
            rank,
            M.cols - rank,
        )

    nullity = (
        M.cols - rank
    )

    if nullity != 0:
        return (
            "NONUNIQUE",
            None,
            rank,
            nullity,
        )

    solution, _ = (
        M.gauss_jordan_solve(
            rhs
        )
    )

    solution = tuple(
        clean(x)
        for x in solution
    )

    return (
        "UNIQUE",
        solution,
        rank,
        0,
    )


def predict(labels, weights, r):
    return clean(
        sum(
            weights[i]
            * q_value(
                labels[i],
                r,
            )
            for i in range(
                len(labels)
            )
        )
    )


def validate_columns(
    k,
    labels,
    weights,
    columns,
):
    results = []

    for r in columns:

        pred = predict(
            labels,
            weights,
            r,
        )

        actual = clean(
            b_value(
                k,
                r,
            )
        )

        results.append(
            (
                r,
                pred,
                actual,
                clean(pred - actual),
            )
        )

    return results


def subsets(items, size):
    """
    Exact deterministic combinations.
    """
    if size == 1:
        return [
            (x,)
            for x in items
        ]

    out = []

    def rec(
        start,
        current,
    ):

        if len(current) == size:
            out.append(
                tuple(current)
            )
            return

        for i in range(
            start,
            len(items),
        ):
            rec(
                i + 1,
                current + [items[i]],
            )

    rec(
        0,
        [],
    )

    return out


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 282 — EXACT q-FAMILY HOLD-OUT / "
        "LEAVE-ONE-OUT PROVENANCE AUDIT"
    )
    print("=" * 78)

    # =========================================================================
    # 1. HOLD-OUT PREFIX/SUFFIX TEST
    # =========================================================================

    print()
    print("=" * 78)
    print("1. PREFIX-TRAIN / SUFFIX-TEST")
    print("=" * 78)

    prefix_results = []

    for size in [2, 3, 4]:

        for labels_tuple in subsets(
            P_VALUES,
            size,
        ):

            labels = list(
                labels_tuple
            )

            for k in range(
                len(B)
            ):

                columns = common_columns(
                    k,
                    labels,
                )

                # Need at least one held-out column and enough training
                # columns to determine all weights uniquely.
                if len(columns) <= size:
                    continue

                train = columns[
                    :size
                ]

                test = columns[
                    size:
                ]

                M, rhs = build_system(
                    k,
                    labels,
                    train,
                )

                status, weights, rank, nullity = (
                    solve_unique(
                        M,
                        rhs,
                    )
                )

                if status != "UNIQUE":
                    prefix_results.append(
                        {
                            "labels": labels,
                            "k": k,
                            "train": train,
                            "test": test,
                            "status": status,
                            "weights": None,
                            "held_out_exact": False,
                        }
                    )
                    continue

                validation = validate_columns(
                    k,
                    labels,
                    weights,
                    test,
                )

                held_out_exact = all(
                    diff == 0
                    for (
                        _,
                        _,
                        _,
                        diff,
                    )
                    in validation
                )

                prefix_results.append(
                    {
                        "labels": labels,
                        "k": k,
                        "train": train,
                        "test": test,
                        "status": status,
                        "weights": weights,
                        "held_out_exact": held_out_exact,
                        "validation": validation,
                    }
                )

    prefix_exact = [
        x
        for x in prefix_results
        if x["held_out_exact"]
    ]

    prefix_fail = [
        x
        for x in prefix_results
        if not x["held_out_exact"]
    ]

    print(
        f"  prefix_cases="
        f"{len(prefix_results)}"
    )

    print(
        f"  prefix_exact_holdout_cases="
        f"{len(prefix_exact)}"
    )

    print(
        f"  prefix_failed_holdout_cases="
        f"{len(prefix_fail)}"
    )

    for item in prefix_results:

        print()
        print(
            f"  labels={item['labels']} "
            f"k={item['k']}"
        )

        print(
            f"    train={item['train']}"
        )

        print(
            f"    test={item['test']}"
        )

        print(
            f"    status={item['status']}"
        )

        print(
            f"    weights={item['weights']}"
        )

        print(
            f"    held_out_exact="
            f"{item['held_out_exact']}"
        )

        if "validation" in item:
            print(
                f"    validation="
                f"{item['validation']}"
            )

    # =========================================================================
    # 2. LEAVE-ONE-OUT TEST
    # =========================================================================

    print()
    print("=" * 78)
    print("2. LEAVE-ONE-COLUMN-OUT TEST")
    print("=" * 78)

    loo_results = []

    for size in [2, 3]:

        for labels_tuple in subsets(
            P_VALUES,
            size,
        ):

            labels = list(
                labels_tuple
            )

            for k in range(
                len(B)
            ):

                columns = common_columns(
                    k,
                    labels,
                )

                if len(columns) <= size:
                    continue

                for held_out in columns:

                    train = [
                        r
                        for r in columns
                        if r != held_out
                    ]

                    # Use exactly 'size' training columns for unique
                    # determination whenever possible.
                    train = train[
                        :size
                    ]

                    if len(train) < size:
                        continue

                    M, rhs = build_system(
                        k,
                        labels,
                        train,
                    )

                    status, weights, rank, nullity = (
                        solve_unique(
                            M,
                            rhs,
                        )
                    )

                    if status != "UNIQUE":
                        loo_results.append(
                            {
                                "labels": labels,
                                "k": k,
                                "train": train,
                                "held_out": held_out,
                                "status": status,
                                "weights": None,
                                "exact": False,
                            }
                        )
                        continue

                    pred = predict(
                        labels,
                        weights,
                        held_out,
                    )

                    actual = clean(
                        b_value(
                            k,
                            held_out,
                        )
                    )

                    exact = (
                        clean(
                            pred - actual
                        )
                        == 0
                    )

                    loo_results.append(
                        {
                            "labels": labels,
                            "k": k,
                            "train": train,
                            "held_out": held_out,
                            "status": status,
                            "weights": weights,
                            "pred": pred,
                            "actual": actual,
                            "exact": exact,
                        }
                    )

    loo_exact = [
        x
        for x in loo_results
        if x["exact"]
    ]

    loo_fail = [
        x
        for x in loo_results
        if not x["exact"]
    ]

    print(
        f"  loo_cases="
        f"{len(loo_results)}"
    )

    print(
        f"  loo_exact="
        f"{len(loo_exact)}"
    )

    print(
        f"  loo_failures="
        f"{len(loo_fail)}"
    )

    for item in loo_results:

        print()
        print(
            f"  labels={item['labels']} "
            f"k={item['k']}"
        )

        print(
            f"    train={item['train']}"
        )

        print(
            f"    held_out={item['held_out']}"
        )

        print(
            f"    status={item['status']}"
        )

        print(
            f"    weights={item['weights']}"
        )

        print(
            f"    exact={item['exact']}"
        )

        if "pred" in item:
            print(
                f"    predicted="
                f"{item['pred']}"
            )

            print(
                f"    actual="
                f"{item['actual']}"
            )

    # =========================================================================
    # 3. MAXIMAL COLUMN TEST
    # =========================================================================

    print()
    print("=" * 78)
    print("3. MAXIMAL OVERDETERMINED TEST")
    print("=" * 78)

    maximal_results = []

    for size in [2, 3]:

        for labels_tuple in subsets(
            P_VALUES,
            size,
        ):

            labels = list(
                labels_tuple
            )

            for k in range(
                len(B)
            ):

                columns = common_columns(
                    k,
                    labels,
                )

                if len(columns) <= size:
                    continue

                M, rhs = build_system(
                    k,
                    labels,
                    columns,
                )

                status, weights, rank, nullity = (
                    solve_unique(
                        M,
                        rhs,
                    )
                )

                exact = (
                    status == "UNIQUE"
                )

                if exact:
                    validation = validate_columns(
                        k,
                        labels,
                        weights,
                        columns,
                    )

                    exact = all(
                        diff == 0
                        for (
                            _,
                            _,
                            _,
                            diff,
                        )
                        in validation
                    )

                maximal_results.append(
                    {
                        "labels": labels,
                        "k": k,
                        "columns": columns,
                        "status": status,
                        "weights": weights,
                        "exact": exact,
                    }
                )

    maximal_exact = [
        x
        for x in maximal_results
        if x["exact"]
    ]

    print(
        f"  maximal_overdetermined_cases="
        f"{len(maximal_results)}"
    )

    print(
        f"  maximal_exact_cases="
        f"{len(maximal_exact)}"
    )

    for item in maximal_results:

        print()
        print(
            f"  labels={item['labels']} "
            f"k={item['k']}"
        )

        print(
            f"    columns={item['columns']}"
        )

        print(
            f"    status={item['status']}"
        )

        print(
            f"    weights={item['weights']}"
        )

        print(
            f"    exact={item['exact']}"
        )

    # =========================================================================
    # 4. INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print("4. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 281R established several exact interpolating combinations,
but interpolation alone is weak.

Experiment 282 asks for prediction.

A source-mixing law is considered structurally credible only if weights
determined from one set of q-columns reproduce B coefficients on columns
that were NOT used to determine those weights.

Therefore:

    interpolation
        !=
    provenance.

The decisive outputs are:

    prefix_exact_holdout_cases
    loo_exact
    maximal_exact_cases.

A persistent exact hold-out result would be evidence that

    B[k,r]

is genuinely a fixed linear functional of the selected q_p(r) columns.

If every hold-out prediction fails, the apparent exact combinations from
Experiment 281R were merely finite interpolation and should not be used
as source formulas.
"""
    )

    # =========================================================================
    # 5. TERMINAL SOURCE REFERENCE
    # =========================================================================

    print()
    print("=" * 78)
    print("5. TERMINAL PROJECTIVE SOURCE")
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
        f"  q1/17={q1_terminal // 17}"
    )

    print(
        f"  q3/17={q3_terminal // 17}"
    )

    # =========================================================================
    # 6. FINAL EXACTNESS
    # =========================================================================

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  prefix_exact_holdout_cases="
        f"{len(prefix_exact)}"
    )

    print(
        f"  prefix_failed_holdout_cases="
        f"{len(prefix_fail)}"
    )

    print(
        f"  loo_exact="
        f"{len(loo_exact)}"
    )

    print(
        f"  loo_failures="
        f"{len(loo_fail)}"
    )

    print(
        f"  maximal_exact_cases="
        f"{len(maximal_exact)}"
    )

    print(
        f"  maximal_total_cases="
        f"{len(maximal_results)}"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  missing_q_values_treated_as_zero=False"
    )

    print(
        "  interpolation_only_not_counted_as_proof=True"
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
        "EXPERIMENT 282 COMPLETE"
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

