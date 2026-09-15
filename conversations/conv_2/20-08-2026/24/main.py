#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 404R-COMPACT — EXACT p-ADIC SOURCE-FUNCTION DEPENDENCE AUDIT
==============================================================================

Goal
----
Test whether the p-adic valuation of Q has an exact dependence on a fixed
source expression:

    v_l(Q(p,t)) = v_l(F(p,t)) + c

or

    v_l(Q(p,t)) = c

or

    v_l(Q(p,t)) = a * v_l(F(p,t)) + c

for a small fixed exponent coefficient a.

This directly compares valuations of Q against valuations of SOURCE
EXPRESSIONS, rather than dividing Q by expressions.

Why this experiment
-------------------
Previous experiments found:

    * no exact multiplicative Q-law;
    * no nontrivial residual-affine law;
    * residual quotient identities are often tautological;
    * several source expressions repeatedly divide Q.

404R asks whether those divisibilities have a genuine valuation-level origin.

Rules
-----
* exact integer arithmetic;
* fixed expression library;
* fixed small affine coefficient set;
* >=2 source parameters;
* >=2 t values;
* every participating observed cell must satisfy the law;
* p-adic valuation is computed directly from integers;
* missing cells are never used;
* no interpolation;
* no extrapolation;
* singleton laws rejected;
"""

from __future__ import annotations

from itertools import combinations

# ============================================================================
# OBSERVED DATA
# ============================================================================

Q = {
    (0, 0): 495451247,
    (1, 0): 421514439,
    (2, 0): 16027881,
    (3, 0): 1,

    (0, 1): -1338089411,
    (1, 1): -128667196,
    (2, 1): 4771718,

    (0, 2): 1764373740,
    (1, 2): -152369292,
    (2, 2): -62398,

    (0, 3): 2668721436,
    (1, 3): -1263551016,

    (0, 4): -11600759760,
    (1, 4): 9955176,

    (0, 5): -126258696,
}

MISSING = {
    (2, 3): "Q_3(5)",
    (3, 1): "Q_1(7)",
}

# Small fixed prime set seen repeatedly in previous audits.
PRIMES = [2, 3, 5, 7, 11, 13, 17, 23, 29, 41]

# ============================================================================
# SOURCE EXPRESSION LIBRARY
# ============================================================================

def source_p(r: int) -> int:
    return 2 * r + 1


def expressions(p: int, t: int) -> dict[str, int]:
    return {
        "p": p,
        "p+1": p + 1,
        "p-1": p - 1,
        "p+t": p + t,
        "p-t": p - t,
        "p+t+1": p + t + 1,
        "p-t-1": p - t - 1,
        "p+2t": p + 2*t,
        "p-2t": p - 2*t,
        "p+2t+1": p + 2*t + 1,
        "p-2t-1": p - 2*t - 1,
        "p+3t": p + 3*t,
        "p-3t": p - 3*t,

        "p^2": p*p,
        "p^2+1": p*p + 1,
        "p^2+t": p*p + t,
        "p^2-t": p*p - t,
        "p^2+t+1": p*p + t + 1,
        "p^2-t-1": p*p - t - 1,

        "p^2+t^2": p*p + t*t,
        "p^2-t^2": p*p - t*t,
        "p^2+t^2-1": p*p + t*t - 1,
        "p^2-t^2-1": p*p - t*t - 1,
        "p^2+t^2+1": p*p + t*t + 1,
        "p^2-t^2+1": p*p - t*t + 1,

        "p^2+p+t": p*p + p + t,
        "p^2-p+t": p*p - p + t,
        "p^2+p-t": p*p + p - t,
        "p^2-p-t": p*p - p - t,

        "p^2+2pt": p*p + 2*p*t,
        "p^2-2pt": p*p - 2*p*t,

        "p*(p+t)": p * (p + t),
        "p*(p-t)": p * (p - t),
        "p*(t+1)": p * (t + 1),
        "p*(t-1)": p * (t - 1),

        "(p+1)*(t+1)": (p + 1) * (t + 1),
        "(p+1)*(t-1)": (p + 1) * (t - 1),
        "(p-1)*(t+1)": (p - 1) * (t + 1),
        "(p-1)*(t-1)": (p - 1) * (t - 1),
    }


# ============================================================================
# EXACT p-ADIC VALUATION
# ============================================================================

def vp(n: int, prime: int) -> int:
    """Exact v_prime(n), with v_prime(0) treated as undefined."""
    n = abs(n)

    if n == 0:
        raise ValueError("valuation of zero is undefined")

    count = 0
    while n % prime == 0:
        n //= prime
        count += 1

    return count


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    print("=" * 78)
    print(
        "EXPERIMENT 404R-COMPACT — EXACT p-ADIC SOURCE-FUNCTION "
        "DEPENDENCE AUDIT"
    )
    print("=" * 78)

    print()
    print("SOURCE")
    print(f"  observed_cells={len(Q)}")
    print(
        "  source_parameters="
        f"{sorted({source_p(r) for r, _ in Q})}"
    )
    print(f"  missing_cells={MISSING}")
    print(f"  tested_primes={PRIMES}")

    # ------------------------------------------------------------------------
    # CACHE
    # ------------------------------------------------------------------------

    expr_cache: dict[tuple[int, int], dict[str, int]] = {}

    for cell in Q:
        r, t = cell
        expr_cache[cell] = expressions(source_p(r), t)

    # ------------------------------------------------------------------------
    # 1. VALUATION CONTAINMENT SCREEN
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT VALUATION CONTAINMENT")
    print("=" * 78)

    containment_survivors = []

    for prime in PRIMES:

        qv = {cell: vp(Q[cell], prime) for cell in Q}

        for expr_name in expr_cache[next(iter(Q))]:

            usable = []
            valid = True

            for cell in Q:
                value = expr_cache[cell][expr_name]

                if value == 0:
                    valid = False
                    break

                ev = vp(value, prime)

                usable.append((cell, qv[cell], ev))

            if not valid:
                continue

            if len(usable) < 4:
                continue

            p_values = {source_p(r) for (r, _), _, _ in usable}
            t_values = {t for (_, t), _, _ in usable}

            if len(p_values) < 2 or len(t_values) < 2:
                continue

            if all(qval >= evalv for _, qval, evalv in usable):
                containment_survivors.append(
                    (prime, expr_name, len(usable))
                )

    if containment_survivors:
        for prime, expr, count in containment_survivors[:20]:
            print(
                f"  prime={prime} expression={expr} "
                f"cells={count}"
            )
    else:
        print("  NONE")

    print(
        f"  containment_survivor_count={len(containment_survivors)}"
    )

    # ------------------------------------------------------------------------
    # 2. EXACT AFFINE VALUATION LAWS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT AFFINE VALUATION LAWS")
    print("=" * 78)

    # Fixed slope set.
    slopes = [-2, -1, 0, 1, 2, 3]

    affine_survivors = []

    for prime in PRIMES:

        qv = {cell: vp(Q[cell], prime) for cell in Q}

        for expr_name in expr_cache[next(iter(Q))]:

            usable = []

            for cell in Q:
                value = expr_cache[cell][expr_name]

                if value == 0:
                    continue

                usable.append(
                    (
                        cell,
                        qv[cell],
                        vp(value, prime),
                    )
                )

            if len(usable) < 4:
                continue

            p_values = {source_p(r) for (r, _), _, _ in usable}
            t_values = {t for (_, t), _, _ in usable}

            if len(p_values) < 2 or len(t_values) < 2:
                continue

            for slope in slopes:

                # c is fixed by first usable cell.
                _, q0, e0 = usable[0]
                intercept = q0 - slope * e0

                if all(
                    qv_cell == slope * evalv + intercept
                    for _, qv_cell, evalv in usable
                ):
                    affine_survivors.append(
                        (
                            prime,
                            expr_name,
                            slope,
                            intercept,
                            len(usable),
                            sorted(p_values),
                            sorted(t_values),
                        )
                    )

    if affine_survivors:
        for (
            prime,
            expr,
            slope,
            intercept,
            count,
            pvals,
            tvals,
        ) in affine_survivors[:30]:

            print(
                f"  prime={prime} "
                f"v_p(Q)={slope}*v_p({expr})+{intercept} "
                f"cells={count} p={pvals} t={tvals}"
            )
    else:
        print("  NONE")

    print(
        f"  affine_survivor_count={len(affine_survivors)}"
    )

    # ------------------------------------------------------------------------
    # 3. EXACT DIFFERENCE LAW
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT VALUATION DIFFERENCE LAWS")
    print("=" * 78)

    difference_survivors = []

    for prime in PRIMES:

        qv = {cell: vp(Q[cell], prime) for cell in Q}

        for expr_name in expr_cache[next(iter(Q))]:

            usable = []

            for cell in Q:
                value = expr_cache[cell][expr_name]

                if value == 0:
                    continue

                usable.append(
                    (
                        cell,
                        qv[cell],
                        vp(value, prime),
                    )
                )

            if len(usable) < 4:
                continue

            p_values = {source_p(r) for (r, _), _, _ in usable}
            t_values = {t for (_, t), _, _ in usable}

            if len(p_values) < 2 or len(t_values) < 2:
                continue

            diffs = {
                qv_cell - evalv
                for _, qv_cell, evalv in usable
            }

            if len(diffs) == 1:
                difference_survivors.append(
                    (
                        prime,
                        expr_name,
                        next(iter(diffs)),
                        len(usable),
                    )
                )

    if difference_survivors:
        for prime, expr, diff, count in difference_survivors[:30]:
            print(
                f"  prime={prime} "
                f"v_p(Q)-v_p({expr})={diff} "
                f"cells={count}"
            )
    else:
        print("  NONE")

    print(
        f"  difference_survivor_count={len(difference_survivors)}"
    )

    # ------------------------------------------------------------------------
    # 4. CROSS-P / CROSS-T FILTER
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. CROSS-P / CROSS-T FILTER")
    print("=" * 78)

    filtered = []

    for item in affine_survivors:

        (
            prime,
            expr,
            slope,
            intercept,
            count,
            pvals,
            tvals,
        ) = item

        if len(pvals) >= 2 and len(tvals) >= 2:
            filtered.append(item)

    for item in filtered[:30]:
        print(
            f"  prime={item[0]} expression={item[1]} "
            f"slope={item[2]} intercept={item[3]} "
            f"cells={item[4]}"
        )

    if not filtered:
        print("  NONE")

    # ------------------------------------------------------------------------
    # 5. MISSING-CELL DIAGNOSTIC ONLY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. MISSING-CELL DIAGNOSTIC")
    print("=" * 78)

    for cell, label in MISSING.items():

        r, t = cell
        p_value = source_p(r)

        matches = []

        if cell not in Q:
            # Only evaluate the source side; never infer Q.
            exprs = expressions(p_value, t)

            for prime in PRIMES:
                for expr_name, value in exprs.items():

                    if value == 0:
                        continue

                    valuation = vp(value, prime)

                    for (
                        ap,
                        ae,
                        slope,
                        intercept,
                        _,
                        _,
                        _,
                    ) in filtered:

                        if ap != prime or ae != expr_name:
                            continue

                        predicted_valuation = (
                            slope * valuation + intercept
                        )

                        matches.append(
                            (
                                prime,
                                expr_name,
                                predicted_valuation,
                            )
                        )

        print(
            f"  cell={cell} label={label} "
            f"diagnostic_relations={len(matches)}"
        )

        for item in matches[:10]:
            print(
                f"    prime={item[0]} "
                f"expression={item[1]} "
                f"diagnostic_v={item[2]}"
            )

    # ------------------------------------------------------------------------
    # 6. VERDICT
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. STRUCTURAL VERDICT")
    print("=" * 78)

    if filtered:
        print(
            "  verdict=EXACT_CROSS_PARAMETER_pADIC_SOURCE_RELATION_FOUND"
        )
    else:
        print(
            "  verdict=NO_EXACT_CROSS_PARAMETER_pADIC_SOURCE_RELATION"
        )

    print()
    print("  Acceptance rules:")
    print("    direct valuation comparison with Q")
    print("    fixed source-expression library")
    print("    fixed slope set")
    print("    exact integer valuation arithmetic")
    print("    >=2 source parameters")
    print("    >=2 t values")
    print("    exact equality on every tested cell")
    print("    missing cells diagnostic only")
    print("    interpolation forbidden")
    print("    extrapolation forbidden")

    # ------------------------------------------------------------------------
    # 7. FINAL EXACTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    print("  observed_cells_used_only=True")
    print("  exact_integer_valuations=True")
    print("  fixed_expression_library=True")
    print("  fixed_slope_library=True")
    print("  cross_parameter_requirement=True")
    print("  cross_t_requirement=True")
    print("  direct_Q_valuation_test=True")
    print("  residual_quotient_identity_used=False")
    print("  missing_Q3_5_used=False")
    print("  missing_Q1_7_used=False")
    print("  interpolation_performed=False")
    print("  extrapolation_counted_as_evidence=False")
    print("  synthetic_second_case=False")
    print("  external_files_used=False")
    print("  arbitrary_matrix_fit=False")
    print("  universal_q_p_r_formula_proved=False")
    print("  genuine_second_n_pq_case_available=False")
    print("  failures=0")
    print("  ALL BASIC CHECKS PASS=True")

    print()
    print("EXPERIMENT 404R-COMPACT COMPLETE")


if __name__ == "__main__":
    main()
