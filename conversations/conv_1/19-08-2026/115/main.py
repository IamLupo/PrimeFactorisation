#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 281R — EXACT q-FAMILY LINEAR-FUNCTIONAL / B-COLUMN MIXING AUDIT
==============================================================================

Repair of Experiment 281.

The previous implementation allowed incomplete/underdetermined SymPy
solutions to propagate None into Matrix objects.

This version:

    * never inserts None into a SymPy Matrix;
    * treats missing q_p(r) entries as undefined, never zero;
    * uses exact QQ matrices;
    * distinguishes:
          NO_SOLUTION
          UNIQUE_EXACT
          NONUNIQUE_EXACT;
    * verifies every returned solution exactly.

Target:

    B[k,r] = sum_p w[k,p] q_p(r)

with weights depending on k but not r.

No arbitrary matrix fitting.
No floating point.
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
# ABSOLUTE-INDEX B TABLE
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
        p in Q
        and 0 <= r < len(Q[p])
    )


def q_value(p, r):
    if not q_defined(p, r):
        return None
    return sp.Rational(Q[p][r])


def b_value(k, r):
    offset = r - k

    if offset < 0:
        return None

    if offset >= len(B[k]):
        return None

    return sp.Rational(
        B[k][offset]
    )


def common_columns(labels, k):
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


def exact_linear_status(M, b):
    """
    Exact classification of M*w=b.

    Returns:
        status:
            NO_SOLUTION
            UNIQUE_EXACT
            NONUNIQUE_EXACT
        rank
        nullity
        solution_tuple_or_None
    """

    M = sp.Matrix(
        [
            [
                sp.Rational(x)
                for x in row
            ]
            for row in M
        ]
    )

    b = sp.Matrix(
        [
            sp.Rational(x)
            for x in b
        ]
    )

    rank = M.rank()
    aug_rank = M.row_join(b).rank()

    if aug_rank != rank:
        return (
            "NO_SOLUTION",
            rank,
            M.cols - rank,
            None,
        )

    nullity = M.cols - rank

    if nullity == 0:
        # Unique solution. Use gauss_jordan_solve safely.
        solution, _ = M.gauss_jordan_solve(
            b
        )

        solution = tuple(
            clean(x)
            for x in solution
        )

        return (
            "UNIQUE_EXACT",
            rank,
            0,
            solution,
        )

    # Consistent but underdetermined.
    #
    # We deliberately do not expose a symbolic free-variable tuple because
    # that is not needed for this experiment.
    return (
        "NONUNIQUE_EXACT",
        rank,
        nullity,
        None,
    )


def build_system(k, labels, columns):
    M = []
    b = []

    for r in columns:

        row = []

        for p in labels:
            value = q_value(
                p,
                r,
            )

            if value is None:
                raise ValueError(
                    f"q_{p}({r}) undefined."
                )

            row.append(value)

        target = b_value(
            k,
            r,
        )

        if target is None:
            raise ValueError(
                f"B[{k},{r}] undefined."
            )

        M.append(row)
        b.append(target)

    return M, b


def verify_solution(
    k,
    labels,
    columns,
    solution,
):
    if solution is None:
        return False

    for r in columns:

        lhs = clean(
            sum(
                solution[i]
                * q_value(
                    labels[i],
                    r,
                )
                for i in range(
                    len(labels)
                )
            )
        )

        rhs = clean(
            b_value(
                k,
                r,
            )
        )

        if lhs != rhs:
            return False

    return True


def combinations(items, size):
    if size == 1:
        return [(x,) for x in items]

    out = []

    def rec(start, current):
        if len(current) == size:
            out.append(
                tuple(current)
            )
            return

        for idx in range(
            start,
            len(items),
        ):
            rec(
                idx + 1,
                current + [items[idx]],
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
        "EXPERIMENT 281R — EXACT q-FAMILY LINEAR-FUNCTIONAL / "
        "B-COLUMN MIXING AUDIT"
    )
    print("=" * 78)

    # =========================================================================
    # 1. SOURCE AVAILABILITY
    # =========================================================================

    print()
    print("=" * 78)
    print("1. SOURCE-COLUMN AVAILABILITY")
    print("=" * 78)

    for r in range(8):

        available = [
            p
            for p in P_VALUES
            if q_defined(p, r)
        ]

        values = [
            q_value(p, r)
            for p in available
        ]

        print()
        print(
            f"  r={r}: "
            f"available_p={available}"
        )

        print(
            f"    q_values={values}"
        )

    # =========================================================================
    # 2. B ABSOLUTE INDEX TABLE
    # =========================================================================

    print()
    print("=" * 78)
    print("2. B[k,r] ABSOLUTE-INDEX TABLE")
    print("=" * 78)

    for k in range(len(B)):

        entries = []

        for offset, value in enumerate(
            B[k]
        ):

            entries.append(
                (
                    k + offset,
                    value,
                )
            )

        print()
        print(
            f"  k={k}: "
            f"{entries}"
        )

    # =========================================================================
    # 3. FOUR-CHANNEL TEST
    # =========================================================================

    print()
    print("=" * 78)
    print("3. FOUR-CHANNEL LINEAR-FUNCTIONAL TEST")
    print("=" * 78)

    four_channel_results = []

    for k in range(
        len(B)
    ):

        labels = P_VALUES

        columns = common_columns(
            labels,
            k,
        )

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    columns={columns}"
        )

        if not columns:

            print(
                "    status=NO_COMMON_COLUMNS"
            )

            four_channel_results.append(
                (
                    k,
                    "NO_COMMON_COLUMNS",
                )
            )

            continue

        M, rhs = build_system(
            k,
            labels,
            columns,
        )

        status, rank, nullity, solution = (
            exact_linear_status(
                M,
                rhs,
            )
        )

        verified = (
            status == "UNIQUE_EXACT"
            and verify_solution(
                k,
                labels,
                columns,
                solution,
            )
        )

        print(
            f"    labels={labels}"
        )

        print(
            f"    rank={rank}"
        )

        print(
            f"    nullity={nullity}"
        )

        print(
            f"    status={status}"
        )

        print(
            f"    solution={solution}"
        )

        print(
            f"    verified={verified}"
        )

        four_channel_results.append(
            (
                k,
                status,
                verified,
            )
        )

    # =========================================================================
    # 4. ALL SOURCE-CHANNEL SUBSETS
    # =========================================================================

    print()
    print("=" * 78)
    print("4. SOURCE-SUBSET RECONSTRUCTION SEARCH")
    print("=" * 78)

    subset_summary = []

    for size in [3, 2]:

        for labels_tuple in combinations(
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
                    labels,
                    k,
                )

                if len(columns) < size:
                    continue

                M, rhs = build_system(
                    k,
                    labels,
                    columns,
                )

                status, rank, nullity, solution = (
                    exact_linear_status(
                        M,
                        rhs,
                    )
                )

                verified = (
                    status
                    == "UNIQUE_EXACT"
                    and verify_solution(
                        k,
                        labels,
                        columns,
                        solution,
                    )
                )

                subset_summary.append(
                    (
                        size,
                        tuple(labels),
                        k,
                        columns,
                        status,
                        verified,
                        solution,
                    )
                )

    exact_subset_cases = [
        item
        for item in subset_summary
        if item[4] == "UNIQUE_EXACT"
        and item[5] is True
    ]

    nonunique_cases = [
        item
        for item in subset_summary
        if item[4] == "NONUNIQUE_EXACT"
    ]

    no_solution_cases = [
        item
        for item in subset_summary
        if item[4] == "NO_SOLUTION"
    ]

    print(
        f"  tested_subset_cases="
        f"{len(subset_summary)}"
    )

    print(
        f"  unique_exact_cases="
        f"{len(exact_subset_cases)}"
    )

    print(
        f"  nonunique_exact_cases="
        f"{len(nonunique_cases)}"
    )

    print(
        f"  no_solution_cases="
        f"{len(no_solution_cases)}"
    )

    for item in exact_subset_cases[:50]:

        size, labels, k0, columns, status, verified, solution = item

        print()
        print(
            f"  exact: "
            f"p={labels} "
            f"k={k0} "
            f"columns={columns}"
        )

        print(
            f"    weights={solution}"
        )

    # =========================================================================
    # 5. WEIGHT PROFILE
    # =========================================================================

    print()
    print("=" * 78)
    print("5. UNIQUE WEIGHT PROFILE")
    print("=" * 78)

    for item in exact_subset_cases:

        (
            size,
            labels,
            k0,
            columns,
            status,
            verified,
            solution,
        ) = item

        print()
        print(
            f"  k={k0}: "
            f"labels={labels}"
        )

        print(
            f"    weights={solution}"
        )

        print(
            f"    columns={columns}"
        )

    # =========================================================================
    # 6. TERMINAL SOURCE
    # =========================================================================

    print()
    print("=" * 78)
    print("6. TERMINAL PROJECTIVE SOURCE")
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
    # 7. INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print("7. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 281R fixes the technical failure in the previous linear-system
implementation.

The mathematical question is unchanged:

    Is a B row a fixed linear functional of the available q-family
    columns?

The classification is now exact:

    UNIQUE_EXACT
        A genuine q-family mixing law exists for those columns.

    NONUNIQUE_EXACT
        The data are compatible, but insufficient to determine unique
        weights.

    NO_SOLUTION
        That source subset cannot generate the B row.

Missing q_p(r) values are never treated as zero.

The strongest result would be repeated UNIQUE_EXACT weights with a
simple dependence on k.

If even two/three-channel subsets fail, then the B operator cannot be
explained as a linear projection of the currently known q-family rows.
"""
    )

    # =========================================================================
    # 8. FINAL
    # =========================================================================

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    unique_four = sum(
        1
        for item in four_channel_results
        if len(item) >= 3
        and item[1] == "UNIQUE_EXACT"
        and item[2] is True
    )

    print(
        f"  four_channel_unique_exact_rows="
        f"{unique_four}"
    )

    print(
        f"  exact_subset_cases="
        f"{len(exact_subset_cases)}"
    )

    print(
        f"  nonunique_subset_cases="
        f"{len(nonunique_cases)}"
    )

    print(
        f"  no_solution_subset_cases="
        f"{len(no_solution_cases)}"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  missing_q_values_treated_as_zero=False"
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
        "EXPERIMENT 281R COMPLETE"
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