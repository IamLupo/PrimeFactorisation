from __future__ import annotations

import sympy as sp


# =============================================================================
# r=5 NEXT EXPERIMENT
# =============================================================================
#
# PURPOSE
# -------
# The previous experiment established exact common support roots.
#
# This experiment moves one level deeper:
#
#       E_j(k,D) = Support_j(D) * Q_j(k,D)
#
# and asks:
#
#   1. At which D-values does the cross-k polynomial degree DROP?
#   2. What are the exact factorizations of those low-degree k-polynomials?
#   3. Do the first nonzero D-slices have a systematic j-dependence?
#   4. Can the k-polynomials be expressed using simple linear factors?
#   5. Can the D-onset and j-dependence be organized into a single
#      boundary table?
#   6. Does a candidate factorial/binomial pattern emerge?
#   7. Are the exceptional j=5 factors genuinely different?
#
# IMPORTANT:
#   - DATA are already exact discrepancies.
#   - No old-law subtraction is performed.
#   - No r=6.
#   - No full pq-kernel expansion.
#   - No replacement universal law is asserted.
#
# The experiment deliberately does NOT trust generic degree-5 interpolation.
# It searches for degree drops and exact factorization at structurally
# distinguished D-values.
# =============================================================================


D, K, X = sp.symbols("D K X")


# =============================================================================
# EXACT SUPPLIED DISCREPANCY DATA
# =============================================================================

DATA: dict[int, dict[int, dict[int, sp.Integer]]] = {
    3: {
        0: {6: -462, 8: -3430, 10: -12684, 12: -50205, 14: -127148, 16: -273294},
        1: {6: 0, 8: -1820, 10: -15315, 12: -44037, 14: -107030, 16: -244790},
        2: {6: 0, 8: -728, 10: -3648, 12: -6930, 14: -53382, 16: -139622},
        3: {6: 0, 8: 0, 10: -525, 12: -14763, 14: -28014, 16: -41706},
        4: {6: 0, 8: 0, 10: -150, 12: 1905, 14: 26814, 16: 0},
        5: {6: 0, 8: 0, 10: 0, 12: 0, 14: -23800, 16: -29892},
    },

    5: {
        0: {6: -2268, 8: -18648, 10: -76125, 12: -219150, 14: -589372, 16: -1272852},
        1: {6: 0, 8: -5712, 10: -34440, 12: -138068, 14: -343728, 16: -755160},
        2: {6: 0, 8: -1632, 10: -10980, 12: -35998, 14: -72888, 16: -291060},
        3: {6: 0, 8: 0, 10: -855, 12: -6255, 14: -53916, 16: -94050},
        4: {6: 0, 8: 0, 10: -190, 12: -1390, 14: 8997, 16: 98560},
        5: {6: 0, 8: 0, 10: 0, 12: 0, 14: 0, 16: -51800},
    },

    7: {
        0: {6: -7524, 8: -60984, 10: -257136, 12: -744996, 14: -1818586, 16: -4072344},
        1: {6: 0, 8: -13860, 10: -83325, 12: -282711, 14: -794430, 16: -1828398},
        2: {6: 0, 8: -3080, 10: -20625, 12: -78298, 14: -158984, 16: -359480},
        3: {6: 0, 8: 0, 10: -1265, 12: -9251, 14: -7546, 16: -66242},
        4: {6: 0, 8: 0, 10: -230, 12: -1682, 14: -735, 16: 135838},
        5: {6: 0, 8: 0, 10: 0, 12: 0, 14: 0, 16: 56056},
    },

    9: {
        0: {6: -19734, 8: -158444, 10: -665236, 12: -1913769, 14: -4553505, 16: -10125050},
        1: {6: 0, 8: -28600, 10: -171600, 12: -570570, 14: -1298440, 16: -3360148},
        2: {6: 0, 8: -5200, 10: -34710, 12: -130221, 14: -189488, 16: -177100},
        3: {6: 0, 8: 0, 10: -1755, 12: -12831, 14: 22750, 16: 533910},
        4: {6: 0, 8: 0, 10: -270, 12: -1974, 14: 3360, 16: 372780},
        5: {6: 0, 8: 0, 10: 0, 12: 0, 14: 0, 16: 132192},
    },

    11: {
        0: {6: -44226, 8: -352716, 10: -1476384, 12: -4138134, 14: -9002994, 16: -20452146},
        1: {6: 0, 8: -52780, 10: -316225, 12: -1028027, 14: -1682590, 16: -3115294},
        2: {6: 0, 8: -8120, 10: -54075, 12: -200417, 14: -12390, 16: 2135770},
        3: {6: 0, 8: 0, 10: -2325, 12: -16995, 14: 96258, 16: 2239830},
        4: {6: 0, 8: 0, 10: -310, 12: -2266, 14: 10672, 16: 1048278},
        5: {6: 0, 8: 0, 10: 0, 12: 0, 14: 0, 16: 277134},
    },

    13: {
        0: {6: -88536, 8: -702576, 10: -2934064, 12: -7980140, 14: -13755924, 16: -27091064},
        1: {6: 0, 8: -89760, 10: -537200, 12: -1702856, 14: -891072, 16: 7361136},
        2: {6: 0, 8: -11968, 10: -79560, 12: -291108, 14: 609960, 16: 11645748},
        3: {6: 0, 8: 0, 10: -2975, 12: -21743, 14: 243984, 16: 6970578},
        4: {6: 0, 8: 0, 10: -350, 12: -2558, 14: 22491, 16: 2485476},
        5: {6: 0, 8: 0, 10: 0, 12: 0, 14: 532252, 16: 533876},
    },
}


KS = sorted(DATA)
JS = list(range(6))
DS = sorted(DATA[KS[0]][0])


# =============================================================================
# BASIC EXACT ACCESS
# =============================================================================

def E(j: int, k: int, d: int) -> sp.Integer:
    return sp.Integer(DATA[k][j][d])


def interpolate_D(j: int, k: int) -> sp.Expr:
    points = [
        (sp.Integer(d), E(j, k, d))
        for d in DS
    ]

    return sp.factor(
        sp.interpolate(points, D)
    )


def common_support_roots(j: int) -> list[int]:
    return [
        d
        for d in DS
        if all(E(j, k, d) == 0 for k in KS)
    ]


def support_factor(j: int) -> sp.Expr:
    s = sp.Integer(1)

    for root in common_support_roots(j):
        s *= D - root

    return sp.factor(s)


def support_quotient_D(j: int, k: int) -> sp.Expr:
    P = sp.Poly(
        interpolate_D(j, k),
        D,
        domain="QQ",
    )

    S = sp.Poly(
        support_factor(j),
        D,
        domain="QQ",
    )

    Q, R = sp.div(P, S)

    if not R.is_zero:
        raise RuntimeError(
            f"support division failed: j={j}, k={k}, "
            f"remainder={sp.factor(R.as_expr())}"
        )

    return sp.factor(Q.as_expr())


def QX(j: int, k: int) -> sp.Expr:
    return sp.factor(
        sp.expand(
            support_quotient_D(j, k).subs(
                D,
                X + 6,
            )
        )
    )


# =============================================================================
# EXACT CROSS-k POLYNOMIAL
# =============================================================================

def cross_k_E(
    j: int,
    d: int,
) -> sp.Expr:

    points = [
        (
            sp.Integer(k),
            E(j, k, d),
        )
        for k in KS
    ]

    return sp.factor(
        sp.interpolate(
            points,
            K,
        )
    )


def cross_k_QX_coefficient(
    j: int,
    power: int,
) -> sp.Expr:

    values = []

    for k in KS:
        p = sp.Poly(
            QX(j, k),
            X,
            domain="QQ",
        )

        values.append(
            (
                sp.Integer(k),
                p.nth(power),
            )
        )

    return sp.factor(
        sp.interpolate(
            values,
            K,
        )
    )


# =============================================================================
# POLYNOMIAL DEGREE HELPERS
# =============================================================================

def polynomial_degree(
    expr: sp.Expr,
    variable,
) -> int:

    expr = sp.expand(expr)

    if expr == 0:
        return -sp.oo

    return int(
        sp.Poly(
            expr,
            variable,
        ).degree()
    )


def factor_over_Q(expr: sp.Expr) -> sp.Expr:
    return sp.factor(
        sp.cancel(
            sp.expand(expr)
        )
    )


# =============================================================================
# FIRST NONZERO D FOR EACH j
# =============================================================================

def first_nonzero_D(
    j: int,
) -> int | None:

    for d in DS:
        if any(
            E(j, k, d) != 0
            for k in KS
        ):
            return d

    return None


# =============================================================================
# ALL CROSS-k DEGREE DROPS
# =============================================================================

def section_1_degree_drop_map() -> None:

    print("=" * 78)
    print("1. CROSS-k DEGREE DROP MAP")
    print("=" * 78)

    for j in JS:

        print(f"j={j}")

        degrees = []

        for d in DS:

            P = cross_k_E(
                j,
                d,
            )

            degree = polynomial_degree(
                P,
                K,
            )

            degrees.append(
                (
                    d,
                    degree,
                )
            )

            print(
                f"  D={d}: "
                f"degree={degree}"
            )

        low = [
            item
            for item in degrees
            if item[1] >= 0 and item[1] < 5
        ]

        print(
            f"  degree-drop slices={low}"
        )

        print()


# =============================================================================
# EXACT FACTORIZATION OF EVERY D-SLICE
# =============================================================================

def section_2_factor_all_D_slices() -> None:

    print("=" * 78)
    print("2. EXACT FACTORIZATION OF CROSS-k D-SLICES")
    print("=" * 78)

    for j in JS:

        print(f"j={j}")

        for d in DS:

            P = cross_k_E(
                j,
                d,
            )

            print(
                f"  D={d}:"
            )
            print(
                f"    degree={polynomial_degree(P,K)}"
            )
            print(
                f"    E(K)={P}"
            )

            if P != 0:
                print(
                    f"    factor(K)="
                    f"{factor_over_Q(P)}"
                )

        print()


# =============================================================================
# STRUCTURALLY DISTINGUISHED D:
#    first nonzero,
#    support roots,
#    first post-support point,
#    terminal point
# =============================================================================

def distinguished_Ds(j: int) -> list[int]:

    result = []

    first = first_nonzero_D(j)

    if first is not None:
        result.append(first)

    result.extend(
        common_support_roots(j)
    )

    roots = common_support_roots(j)

    if roots:
        after = [
            d
            for d in DS
            if d > max(roots)
        ]

        if after:
            result.append(after[0])

    if DS:
        result.append(DS[-1])

    return sorted(
        set(result)
    )


# =============================================================================
# DISTINGUISHED-SLICE FACTORIZATION
# =============================================================================

def section_3_distinguished_slices() -> None:

    print("=" * 78)
    print("3. DISTINGUISHED D-SLICE FACTORIZATION")
    print("=" * 78)

    for j in JS:

        print(f"j={j}")
        print(
            f"  support={support_factor(j)}"
        )
        print(
            f"  roots={common_support_roots(j)}"
        )
        print(
            f"  first_nonzero_D="
            f"{first_nonzero_D(j)}"
        )

        for d in distinguished_Ds(j):

            P = cross_k_E(
                j,
                d,
            )

            print(
                f"  D={d}:"
            )
            print(
                f"    degree={polynomial_degree(P,K)}"
            )
            print(
                f"    factor={factor_over_Q(P)}"
            )

        print()


# =============================================================================
# LINEAR FACTOR AUDIT
# =============================================================================

def linear_factor_roots(
    expr: sp.Expr,
) -> list[sp.Expr]:

    if expr == 0:
        return []

    P = sp.Poly(
        expr,
        K,
        domain="QQ",
    )

    roots = sp.roots(
        P.as_expr(),
        K,
    )

    return sorted(
        list(roots.keys()),
        key=lambda x: str(x),
    )


def section_4_k_root_audit() -> None:

    print("=" * 78)
    print("4. EXACT k-ROOT AUDIT FOR DEGREE-DROP SLICES")
    print("=" * 78)

    for j in JS:

        print(f"j={j}")

        for d in DS:

            P = cross_k_E(
                j,
                d,
            )

            degree = polynomial_degree(
                P,
                K,
            )

            if degree == -sp.oo:
                continue

            if degree < 5:

                print(
                    f"  D={d}:"
                )
                print(
                    f"    degree={degree}"
                )
                print(
                    f"    factor="
                    f"{factor_over_Q(P)}"
                )
                print(
                    f"    exact roots="
                    f"{linear_factor_roots(P)}"
                )

        print()


# =============================================================================
# FIRST-NONZERO BOUNDARY TABLE
# =============================================================================

def section_5_boundary_table() -> None:

    print("=" * 78)
    print("5. FIRST-NONZERO D BOUNDARY TABLE")
    print("=" * 78)

    for j in JS:

        d0 = first_nonzero_D(j)

        if d0 is None:
            print(
                f"j={j}: no nonzero D in grid"
            )
            continue

        P = cross_k_E(
            j,
            d0,
        )

        print(
            f"j={j}: "
            f"D0={d0}"
        )
        print(
            f"  degree={polynomial_degree(P,K)}"
        )
        print(
            f"  E(K)={P}"
        )
        print(
            f"  factor={factor_over_Q(P)}"
        )
        print(
            f"  roots={linear_factor_roots(P)}"
        )
        print()


# =============================================================================
# FIRST NONZERO D COMBINATORIAL PATTERN SEARCH
# =============================================================================

def section_6_boundary_factor_patterns() -> None:

    print("=" * 78)
    print("6. BOUNDARY FACTOR PATTERN SEARCH")
    print("=" * 78)

    for j in JS:

        d0 = first_nonzero_D(j)

        if d0 is None:
            continue

        P = cross_k_E(
            j,
            d0,
        )

        print(f"j={j}, D0={d0}")
        print(
            f"  factor(K)="
            f"{factor_over_Q(P)}"
        )

        candidates = {
            "K+1": K + 1,
            "K+2": K + 2,
            "K+3": K + 3,
            "K+4": K + 4,
            "K+5": K + 5,
            "2K+5": 2*K + 5,
            "2K+7": 2*K + 7,
            "2K+9": 2*K + 9,
            "2K+11": 2*K + 11,
        }

        print("  candidate exact divisibility:")

        for name, factor in candidates.items():

            q, r = sp.div(
                sp.Poly(P, K),
                sp.Poly(factor, K),
            )

            if r.is_zero:
                print(
                    f"    {name}: YES"
                )
            else:
                print(
                    f"    {name}: no"
                )

        print()


# =============================================================================
# SEARCH PRODUCTS OF CONSECUTIVE LINEAR FACTORS
# =============================================================================

def factor_candidate_products(
    max_shift: int = 8,
) -> list[tuple[str, sp.Expr]]:

    candidates = []

    for length in range(1, 7):

        for start in range(
            -2,
            max_shift + 1,
        ):

            factors = [
                K + start + i
                for i in range(length)
            ]

            expr = sp.prod(
                factors
            )

            candidates.append(
                (
                    "*".join(
                        str(f)
                        for f in factors
                    ),
                    expr,
                )
            )

    return candidates


def section_7_consecutive_factor_search() -> None:

    print("=" * 78)
    print("7. CONSECUTIVE-LINEAR-FACTOR SEARCH")
    print("=" * 78)

    candidates = factor_candidate_products()

    extra = [
        ("2K+5", 2*K + 5),
        ("2K+7", 2*K + 7),
        ("2K+9", 2*K + 9),
        ("2K+11", 2*K + 11),
        ("2K+13", 2*K + 13),
    ]

    candidates.extend(extra)

    for j in JS:

        d0 = first_nonzero_D(j)

        if d0 is None:
            continue

        P = sp.Poly(
            cross_k_E(j, d0),
            K,
            domain="QQ",
        )

        print(
            f"j={j}, D0={d0}"
        )
        print(
            f"  P(K)={factor_over_Q(P.as_expr())}"
        )

        hits = []

        for name, candidate in candidates:

            C = sp.Poly(
                candidate,
                K,
                domain="QQ",
            )

            q, r = sp.div(
                P,
                C,
            )

            if r.is_zero:
                hits.append(
                    (
                        name,
                        sp.factor(q.as_expr()),
                    )
                )

        if hits:
            for name, quotient in hits:
                print(
                    f"  divides by {name}: "
                    f"quotient={quotient}"
                )
        else:
            print(
                "  no tested consecutive-factor product divides"
            )

        print()


# =============================================================================
# SUPPORT-REMOVED FIRST-POST-ROOT SLICE
# =============================================================================

def section_8_first_post_support() -> None:

    print("=" * 78)
    print("8. FIRST POST-SUPPORT D-SLICE")
    print("=" * 78)

    for j in JS:

        roots = common_support_roots(j)

        if not roots:
            continue

        post = [
            d
            for d in DS
            if d > max(roots)
        ]

        if not post:
            continue

        d = post[0]

        P = cross_k_E(
            j,
            d,
        )

        print(
            f"j={j}: "
            f"post-support D={d}"
        )
        print(
            f"  degree={polynomial_degree(P,K)}"
        )
        print(
            f"  factor={factor_over_Q(P)}"
        )
        print(
            f"  roots={linear_factor_roots(P)}"
        )
        print()


# =============================================================================
# TERMINAL D=16
# =============================================================================

def section_9_terminal_D16() -> None:

    print("=" * 78)
    print("9. TERMINAL D=16 CROSS-k STRUCTURE")
    print("=" * 78)

    for j in JS:

        P = cross_k_E(
            j,
            16,
        )

        print(
            f"j={j}: "
            f"degree={polynomial_degree(P,K)}"
        )
        print(
            f"  E(K)={P}"
        )
        print(
            f"  factor={factor_over_Q(P)}"
        )
        print(
            f"  roots={linear_factor_roots(P)}"
        )
        print()


# =============================================================================
# J=5 D=14 AND D=16 SIDE BY SIDE
# =============================================================================

def section_10_j5_terminal_pair() -> None:

    print("=" * 78)
    print("10. J=5 D=14 / D=16 SIDE-BY-SIDE FACTORIZATION")
    print("=" * 78)

    j = 5

    for d in [14, 16]:

        P = cross_k_E(
            j,
            d,
        )

        print(
            f"D={d}"
        )
        print(
            f"  degree={polynomial_degree(P,K)}"
        )
        print(
            f"  E(K)={P}"
        )
        print(
            f"  factor={factor_over_Q(P)}"
        )
        print(
            f"  roots={linear_factor_roots(P)}"
        )
        print()


# =============================================================================
# COMPARE SUCCESSIVE j VALUES AT SAME D
# =============================================================================

def section_11_j_transition() -> None:

    print("=" * 78)
    print("11. SUCCESSIVE-j TRANSITION TABLE")
    print("=" * 78)

    for d in DS:

        print(
            f"D={d}"
        )

        previous = None

        for j in JS:

            P = cross_k_E(
                j,
                d,
            )

            factored = factor_over_Q(
                P
            )

            print(
                f"  j={j}: "
                f"degree={polynomial_degree(P,K)} "
                f"factor={factored}"
            )

            if previous is not None:

                ratio = None

                if previous != 0 and P != 0:
                    ratio = sp.factor(
                        sp.cancel(
                            P / previous
                        )
                    )

                print(
                    f"    ratio_to_previous_j="
                    f"{ratio}"
                )

            previous = P

        print()


# =============================================================================
# FIXED-D FINITE DIFFERENCES IN k
# =============================================================================

def finite_difference_degree(
    values: list[sp.Expr],
) -> int:

    row = [
        sp.factor(v)
        for v in values
    ]

    if all(v == 0 for v in row):
        return -sp.oo

    for degree in range(
        len(row)
    ):

        if len(row) == 1:
            return degree

        if len(set(row)) == 1:
            if row[0] == 0:
                return degree - 1

        row = [
            sp.factor(
                row[i + 1] - row[i]
            )
            for i in range(
                len(row) - 1
            )
        ]

    return len(values) - 1


def section_12_fixed_D_difference_degrees() -> None:

    print("=" * 78)
    print("12. EXACT FINITE-DIFFERENCE DEGREE IN k AT EACH D")
    print("=" * 78)

    for j in JS:

        print(
            f"j={j}"
        )

        for d in DS:

            values = [
                E(j, k, d)
                for k in KS
            ]

            degree = finite_difference_degree(
                values
            )

            print(
                f"  D={d}: "
                f"degree={degree}"
            )

        print()


# =============================================================================
# FACTOR THE FIRST NONZERO SLICE IN TERMS OF SHIFTED K
# =============================================================================

def section_13_shifted_K_boundary() -> None:

    print("=" * 78)
    print("13. SHIFTED-K FACTORIZATION OF FIRST-NONZERO SLICES")
    print("=" * 78)

    for j in JS:

        d0 = first_nonzero_D(j)

        if d0 is None:
            continue

        P = cross_k_E(
            j,
            d0,
        )

        shifts = [
            ("K", K),
            ("K+1", K+1),
            ("K+2", K+2),
            ("K+3", K+3),
            ("K+4", K+4),
            ("K+5", K+5),
        ]

        print(
            f"j={j}, D0={d0}"
        )

        for name, Y in shifts:

            # K = Y - shift
            if name == "K":
                shift = 0
            elif name == "K+1":
                shift = 1
            elif name == "K+2":
                shift = 2
            elif name == "K+3":
                shift = 3
            elif name == "K+4":
                shift = 4
            else:
                shift = 5

            shifted = sp.factor(
                sp.expand(
                    P.subs(
                        K,
                        X - shift,
                    )
                )
            )

            print(
                f"  variable {name}: "
                f"{shifted}"
            )

        print()


# =============================================================================
# SEARCH FOR A GENERAL FIRST-NONZERO FORM
#
# Try:
#
#   C * product(K+a_i) * (2K+b)
#
# with small integer shifts.
# =============================================================================

def section_14_boundary_template_search() -> None:

    print("=" * 78)
    print("14. FIRST-NONZERO TEMPLATE SEARCH")
    print("=" * 78)

    # Candidate product lengths.
    for j in JS:

        d0 = first_nonzero_D(j)

        if d0 is None:
            continue

        P = sp.Poly(
            cross_k_E(j, d0),
            K,
            domain="QQ",
        )

        degree = P.degree()

        print(
            f"j={j}, D0={d0}, "
            f"degree={degree}"
        )

        if degree < 1:
            print()
            continue

        expr = P.as_expr()

        # Rational-root / exact linear-factor decomposition.
        roots = sp.solve(
            sp.Eq(
                expr,
                0,
            ),
            K,
        )

        print(
            f"  exact rational/algebraic roots="
            f"{roots}"
        )

        # Factor integer-shift factors K+a.
        detected = []

        for a in range(
            -15,
            16,
        ):

            factor = K + a

            q, r = sp.div(
                P,
                sp.Poly(
                    factor,
                    K,
                    domain="QQ",
                ),
            )

            count = 0
            current = P

            while True:

                quotient, remainder = sp.div(
                    current,
                    sp.Poly(
                        factor,
                        K,
                        domain="QQ",
                    ),
                )

                if not remainder.is_zero:
                    break

                count += 1
                current = quotient

            if count > 0:
                detected.append(
                    (
                        f"K+({a})",
                        count,
                    )
                )

        print(
            f"  detected integer-shift factors="
            f"{detected}"
        )

        print()


# =============================================================================
# SUPPORT AND FIRST-NONZERO D COMPARISON
# =============================================================================

def section_15_support_vs_onset() -> None:

    print("=" * 78)
    print("15. SUPPORT ROOTS VS FIRST-NONZERO D")
    print("=" * 78)

    for j in JS:

        roots = common_support_roots(j)
        d0 = first_nonzero_D(j)

        print(
            f"j={j}: "
            f"support_roots={roots}, "
            f"first_nonzero_D={d0}"
        )

    print()


# =============================================================================
# OVERLAP OF FACTORS BETWEEN j VALUES
# =============================================================================

def section_16_factor_overlap() -> None:

    print("=" * 78)
    print("16. FACTOR OVERLAP BETWEEN SUCCESSIVE j")
    print("=" * 78)

    for d in DS:

        print(
            f"D={d}"
        )

        for j1 in JS:
            for j2 in JS:
                if j2 <= j1:
                    continue

                P1 = sp.Poly(
                    cross_k_E(j1, d),
                    K,
                    domain="QQ",
                )

                P2 = sp.Poly(
                    cross_k_E(j2, d),
                    K,
                    domain="QQ",
                )

                if P1.is_zero or P2.is_zero:
                    gcd = 0
                else:
                    gcd = sp.factor(
                        sp.gcd(
                            P1,
                            P2,
                        ).as_expr()
                    )

                if gcd not in [1, -1, 0]:
                    print(
                        f"  j=({j1},{j2}): "
                        f"gcd={gcd}"
                    )

        print()


# =============================================================================
# RECONSTRUCTION CHECK
# =============================================================================

def section_17_reconstruction() -> None:

    print("=" * 78)
    print("17. EXACT RECONSTRUCTION CHECK")
    print("=" * 78)

    tested = 0
    failures = 0

    for j in JS:

        for k in KS:

            P = interpolate_D(
                j,
                k,
            )

            for d in DS:

                tested += 1

                if sp.factor(
                    P.subs(
                        D,
                        d,
                    )
                    -
                    E(j, k, d)
                ) != 0:

                    failures += 1

                    print(
                        "FAIL:",
                        j,
                        k,
                        d,
                    )

    print(
        f"tested={tested}"
    )
    print(
        f"failures={failures}"
    )
    print()


# =============================================================================
# KEY STRUCTURAL TABLE
# =============================================================================

def section_18_key_table() -> None:

    print("=" * 78)
    print("18. KEY STRUCTURAL TABLE")
    print("=" * 78)

    print(
        "j | support | support_deg | first_nonzero_D | "
        "first_slice_degree | first_slice_factor"
    )
    print(
        "-" * 78
    )

    for j in JS:

        support = support_factor(j)
        d0 = first_nonzero_D(j)

        if d0 is None:
            degree = None
            factor = None
        else:
            P = cross_k_E(
                j,
                d0,
            )

            degree = polynomial_degree(
                P,
                K,
            )

            factor = factor_over_Q(
                P
            )

        print(
            f"{j} | "
            f"{support} | "
            f"{sp.Poly(support,D).degree()} | "
            f"{d0} | "
            f"{degree} | "
            f"{factor}"
        )

    print()


# =============================================================================
# FINAL DIAGNOSTIC
# =============================================================================

def section_19_final_diagnostic() -> None:

    print("=" * 78)
    print("19. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "This experiment deliberately avoids fitting a new universal law."
    )
    print()

    print(
        "The important object is the exact fixed-D cross-k polynomial:"
    )
    print()
    print(
        "    P_{j,D}(K) = E_j(K,D)"
    )
    print()

    print(
        "A degree below 5 at a specific D is potentially meaningful "
        "because six k-values are available."
    )
    print()

    print(
        "The first nonzero D slice is especially important because it "
        "isolates the onset of the correction."
    )
    print()

    print(
        "The script therefore records:"
    )
    print(
        "  - exact degree in k;"
    )
    print(
        "  - exact rational factorization;"
    )
    print(
        "  - exact k-roots;"
    )
    print(
        "  - repeated integer-shift factors;"
    )
    print(
        "  - overlap of factors between j-values;"
    )
    print(
        "  - first-nonzero boundary behaviour."
    )
    print()

    print(
        "The previous output already shows striking low-degree onset slices:"
    )
    print(
        "  j=0, D=6"
    )
    print(
        "  j=1, D=8"
    )
    print(
        "  j=2, D=8"
    )
    print(
        "  j=3, D=10"
    )
    print(
        "  j=4, D=10"
    )
    print(
        "  j=5, D=14"
    )
    print()

    print(
        "Do NOT infer a law from those six values alone."
    )
    print(
        "Use this experiment to identify exact factor structure that "
        "can subsequently be tested on fresh k,D data."
    )
    print()

    print(
        "No r=6."
    )
    print(
        "No full pq-kernel expansion."
    )
    print(
        "No replacement universal r,j formula."
    )
    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    print("=" * 78)
    print("R=5 NEXT EXPERIMENT: CROSS-k BOUNDARY FACTORIZATION")
    print("=" * 78)
    print()

    print(
        f"k values = {KS}"
    )
    print(
        f"j values = {JS}"
    )
    print(
        f"D values = {DS}"
    )
    print()

    section_1_degree_drop_map()
    section_2_factor_all_D_slices()
    section_3_distinguished_slices()
    section_4_k_root_audit()
    section_5_boundary_table()
    section_6_boundary_factor_patterns()
    section_7_consecutive_factor_search()
    section_8_first_post_support()
    section_9_terminal_D16()
    section_10_j5_terminal_pair()
    section_11_j_transition()
    section_12_fixed_D_difference_degrees()
    section_13_shifted_K_boundary()
    section_14_boundary_template_search()
    section_15_support_vs_onset()
    section_16_factor_overlap()
    section_17_reconstruction()
    section_18_key_table()
    section_19_final_diagnostic()


if __name__ == "__main__":
    main()
