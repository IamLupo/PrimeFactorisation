#!/usr/bin/env python3

"""
==============================================================================
CYCLOTOMIC RESULTANT / GCD / ORDER-STRUCTURE EXPERIMENT
==============================================================================

Purpose
-------

The previous experiment established:

    Phi_3(r) = Phi_6(r+1)

and showed that blind searching for small roots does not scale.

The next experiment therefore studies the algebraic interaction between

    Phi_3(x) = x^2 + x + 1
    Phi_6(x) = x^2 - x + 1

and the n-only quantities

    Phi_3(n) = n^2+n+1
    Phi_6(n) = n^2-n+1.

The main questions are:

1. What are the exact gcd/resultant identities between
       Phi_i(r)
   and
       Phi_j(n)?

2. What happens after substitutions such as

       r+n
       r-n
       r*n
       r*(r+1)
       r*(r-1)
       n+r+1
       n-r-1

3. Are there identities forcing a common divisor of
       Phi_i(r)
   and
       Phi_j(n)
   to divide a simple expression in n or r?

4. Can the observed gcds be explained completely by cyclotomic
   identities, or is there hidden factorization information?

5. Does the multiplicative-order interpretation create an n-only
   constraint that can predict a useful auxiliary r?

6. Can we derive a polynomial whose roots encode the hidden
   cubic roots modulo p or q?

7. Can gcds involving n, Phi_3(r), Phi_6(r), Phi_3(n), Phi_6(n)
   produce an actual nontrivial divisor of n?

NO CSV FILES.
EVERYTHING IS PRINTED.

==============================================================================
"""

from __future__ import annotations

from math import gcd
import sympy as sp


# =============================================================================
# CONFIGURATION
# =============================================================================

REP_P = 22612043
REP_Q = 26706517
REP_N = REP_P * REP_Q

SMALL_RS = list(range(1, 41))

SMALL_MODULI = [
    3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47
]

# Keep symbolic searches manageable.
R_SYMBOLS = list(range(1, 16))

K_VALUES = [3, 6]


# =============================================================================
# BASIC POLYNOMIALS
# =============================================================================

def phi3(x):
    return x**2 + x + 1


def phi6(x):
    return x**2 - x + 1


def bitlen(x):
    return abs(int(x)).bit_length() if x else 1


def section(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def line():
    print("-" * 78)


# =============================================================================
# SYMBOLIC VARIABLES
# =============================================================================

r, s, n, u, x = sp.symbols(
    "r s n u x"
)

P3 = phi3(r)
P6 = phi6(r)

N3 = phi3(n)
N6 = phi6(n)


# =============================================================================
# SECTION 1
# EXACT SYMBOLIC RELATIONS
# =============================================================================

def test_exact_identities():

    section(
        "1. EXACT CYCLOTOMIC IDENTITIES"
    )

    print()
    print(
        "Phi_3(r) =",
        sp.expand(phi3(r))
    )

    print(
        "Phi_6(r) =",
        sp.expand(phi6(r))
    )

    print()

    identities = {
        "Phi3(r)-Phi6(r+1)":
            phi3(r) - phi6(r + 1),

        "Phi3(r)-Phi6(r-1)":
            phi3(r) - phi6(r - 1),

        "Phi6(r)-Phi3(r-1)":
            phi6(r) - phi3(r - 1),

        "Phi6(r)-Phi3(r+1)":
            phi6(r) - phi3(r + 1),
    }

    for name, expr in identities.items():

        print(
            f"{name:32s} = {sp.factor(expr)}"
        )

    print()

    print(
        "Expected exact relation:"
    )

    print(
        "Phi_3(r) = Phi_6(r+1)"
    )


# =============================================================================
# SECTION 2
# RESULTANTS
# =============================================================================

def test_resultants():

    section(
        "2. RESULTANTS BETWEEN CYCLOTOMIC FAMILIES"
    )

    print()
    print(
        "We eliminate r and compute exact resultants."
    )

    polys = {
        "Phi3(r)": phi3(r),
        "Phi6(r)": phi6(r),
        "Phi3(r+1)": phi3(r + 1),
        "Phi6(r+1)": phi6(r + 1),
        "Phi3(r-1)": phi3(r - 1),
        "Phi6(r-1)": phi6(r - 1),
    }

    names = list(polys.keys())

    print()
    print(
        "Resultant(f,g) with respect to r"
    )

    line()

    for i in range(len(names)):

        for j in range(i + 1, len(names)):

            a = names[i]
            b = names[j]

            result = sp.factor(
                sp.resultant(
                    polys[a],
                    polys[b],
                    r
                )
            )

            print(
                f"Res({a:12s}, {b:12s}) = {result}"
            )


# =============================================================================
# SECTION 3
# GCD OF SYMBOLIC POLYNOMIALS
# =============================================================================

def test_symbolic_gcds():

    section(
        "3. SYMBOLIC GCD / EUCLIDEAN RELATIONS"
    )

    expressions = {
        "Phi3(r)": phi3(r),
        "Phi6(r)": phi6(r),
        "Phi3(n)": phi3(n),
        "Phi6(n)": phi6(n),
    }

    print()

    for a_name, a in expressions.items():

        for b_name, b in expressions.items():

            if a_name >= b_name:
                continue

            print(
                f"gcd({a_name}, {b_name}) =",
                sp.factor(
                    sp.gcd(
                        sp.Poly(a, r, n),
                        sp.Poly(b, r, n)
                    ).as_expr()
                )
            )


# =============================================================================
# SECTION 4
# GCD PATTERNS BETWEEN Phi(r) AND Phi(n)
# =============================================================================

def test_numeric_gcd_table():

    section(
        "4. NUMERIC GCD TABLE: Phi_i(r) vs Phi_j(n)"
    )

    N3 = phi3(REP_N)
    N6 = phi6(REP_N)

    print()
    print(
        f"n = {REP_N}"
    )

    print()
    print(
        "r      Phi3(r) gcd Phi3(n)"
        "     Phi3(r) gcd Phi6(n)"
        "     Phi6(r) gcd Phi3(n)"
        "     Phi6(r) gcd Phi6(n)"
    )

    line()

    for rv in SMALL_RS:

        a = phi3(rv)
        b = phi6(rv)

        g33 = gcd(a, N3)
        g36 = gcd(a, N6)
        g63 = gcd(b, N3)
        g66 = gcd(b, N6)

        print(
            f"{rv:2d}"
            f"{g33:24d}"
            f"{g36:24d}"
            f"{g63:24d}"
            f"{g66:24d}"
        )


# =============================================================================
# SECTION 5
# FACTOR THE GCD PATTERNS
# =============================================================================

def factor_gcd_patterns():

    section(
        "5. FACTORIZATION OF NONTRIVIAL GCD PATTERNS"
    )

    N3 = phi3(REP_N)
    N6 = phi6(REP_N)

    print()
    print(
        "Only nontrivial gcds are shown."
    )

    line()

    for rv in SMALL_RS:

        values = {
            "gcd(Phi3(r),Phi3(n))":
                gcd(phi3(rv), N3),

            "gcd(Phi3(r),Phi6(n))":
                gcd(phi3(rv), N6),

            "gcd(Phi6(r),Phi3(n))":
                gcd(phi6(rv), N3),

            "gcd(Phi6(r),Phi6(n))":
                gcd(phi6(rv), N6),
        }

        for name, g in values.items():

            if g > 1:

                print(
                    f"r={rv:2d}"
                    f" {name:35s}"
                    f" = {g}"
                    f" = {sp.factorint(g)}"
                )


# =============================================================================
# SECTION 6
# SUBSTITUTION SEARCH
# =============================================================================

def substitution_search():

    section(
        "6. SUBSTITUTION SEARCH"
    )

    substitutions = {
        "r+n": r + n,
        "r-n": r - n,
        "n+r+1": n + r + 1,
        "n-r-1": n - r - 1,
        "r*(r+1)": r * (r + 1),
        "r*(r-1)": r * (r - 1),
        "r^2": r**2,
        "n*r": n * r,
        "n+r": n + r,
        "n-r": n - r,
    }

    base = {
        "Phi3(r)": phi3(r),
        "Phi6(r)": phi6(r),
    }

    n_polys = {
        "Phi3(n)": phi3(n),
        "Phi6(n)": phi6(n),
    }

    print()
    print(
        "Search expressions gcd(Phi_i(substitution), Phi_j(n))."
    )

    for sub_name, sub in substitutions.items():

        print()
        print(
            f"SUBSTITUTION: {sub_name}"
        )

        for phi_name, base_poly in base.items():

            transformed = (
                base_poly.subs(
                    r,
                    sub
                )
            )

            for n_name, n_poly in n_polys.items():

                try:

                    result = sp.factor(
                        sp.resultant(
                            transformed,
                            n_poly,
                            n
                        )
                    )

                    print(
                        f"Res({phi_name}({sub_name}),"
                        f" {n_name}) = {result}"
                    )

                except Exception as exc:

                    print(
                        f"Could not compute:"
                        f" {phi_name}({sub_name}) vs {n_name}"
                    )

                    print(
                        f"Reason: {exc}"
                    )


# =============================================================================
# SECTION 7
# GCDS AFTER SIMPLE TRANSFORMATIONS
# =============================================================================

def numerical_substitution_gcds():

    section(
        "7. NUMERICAL TRANSFORMATION GCD SEARCH"
    )

    N3 = phi3(REP_N)
    N6 = phi6(REP_N)

    substitutions = {
        "r+n":
            lambda rv: rv + REP_N,

        "r-n":
            lambda rv: rv - REP_N,

        "n+r+1":
            lambda rv: REP_N + rv + 1,

        "n-r-1":
            lambda rv: REP_N - rv - 1,

        "r*(r+1)":
            lambda rv: rv * (rv + 1),

        "r*(r-1)":
            lambda rv: rv * (rv - 1),

        "r^2":
            lambda rv: rv * rv,
    }

    for name, fn in substitutions.items():

        hits = []

        for rv in SMALL_RS:

            z = fn(rv)

            g1 = gcd(
                phi3(z),
                N3
            )

            g2 = gcd(
                phi3(z),
                N6
            )

            g3 = gcd(
                phi6(z),
                N3
            )

            g4 = gcd(
                phi6(z),
                N6
            )

            if max(
                g1, g2, g3, g4
            ) > 1:

                hits.append(
                    (
                        rv,
                        g1,
                        g2,
                        g3,
                        g4
                    )
                )

        print()
        print(
            f"Transformation: {name}"
        )

        if not hits:

            print(
                "  no nontrivial hits"
            )

        else:

            for item in hits:

                print(
                    f"  r={item[0]}"
                    f" gcds={item[1:]}"
                )


# =============================================================================
# SECTION 8
# ORDER INTERPRETATION
# =============================================================================

def order_interpretation():

    section(
        "8. MULTIPLICATIVE-ORDER INTERPRETATION"
    )

    print()
    print(
        "Phi3(r)=0 mod p means r has order 3 modulo p."
    )

    print(
        "Phi6(r)=0 mod p means r has order 6 modulo p."
    )

    print()
    print(
        "For prime p != 3, a nontrivial cubic root exists iff"
    )

    print(
        "p == 1 mod 3."
    )

    print()
    print(
        "Representative factorization:"
    )

    print(
        f"p={REP_P}, p mod 3={REP_P % 3}"
    )

    print(
        f"q={REP_Q}, q mod 3={REP_Q % 3}"
    )

    print()

    for factor, name in [
        (REP_P, "p"),
        (REP_Q, "q"),
    ]:

        print(
            f"{name}:"
        )

        print(
            f"  factor mod 3 = {factor % 3}"
        )

        roots3 = [
            rr
            for rr in range(
                min(factor, 1000)
            )
            if phi3(rr) % factor == 0
        ]

        # For large factors we use direct modular solving.
        if factor > 1000:

            roots3 = sp.sqrt_mod(
                (-3) % factor,
                factor,
                all_roots=True
            )

            roots3 = [
                (
                    -1 + z
                ) * pow(
                    2,
                    -1,
                    factor
                ) % factor
                for z in roots3
            ]

            roots3 = sorted(
                set(roots3)
            )

        print(
            f"  Phi3 roots = {roots3}"
        )

        print()


# =============================================================================
# SECTION 9
# HIDDEN-FACTOR GCD SEARCH USING RESULTANT CANDIDATES
# =============================================================================

def resultant_candidate_search():

    section(
        "9. RESULTANT-INSPIRED FACTOR SEARCH"
    )

    print()
    print(
        "For each small r, construct candidate integers from"
    )

    print(
        "gcd/Phi/resultant-style expressions and test whether"
    )

    print(
        "they produce a divisor of n."
    )

    print()
    print(
        "Only n and r are used."
    )

    candidates = []

    for rv in SMALL_RS:

        a = phi3(rv)
        b = phi6(rv)

        vals = {
            "Phi3(r)": a,
            "Phi6(r)": b,
            "Phi3(r)*Phi6(r)": a * b,
            "Phi3(r)-Phi6(r)": a - b,
            "Phi3(r)+Phi6(r)": a + b,
            "Phi3(r)*Phi6(r)-1":
                a * b - 1,
            "Phi3(r)*Phi6(r)+1":
                a * b + 1,
        }

        for name, value in vals.items():

            g = gcd(
                REP_N,
                abs(value)
            )

            if 1 < g < REP_N:

                candidates.append(
                    (
                        rv,
                        name,
                        value,
                        g
                    )
                )

    if not candidates:

        print(
            "No direct n-factor found."
        )

    else:

        for rv, name, value, g in candidates:

            print(
                f"r={rv}"
                f" {name}"
                f" value={value}"
                f" gcd={g}"
            )


# =============================================================================
# SECTION 10
# DUAL PHI RESULTANTS WITH FIXED n
# =============================================================================

def dual_n_resultants():

    section(
        "10. DUAL RESULTANT SEARCH WITH n FIXED"
    )

    print()
    print(
        "Compute resultants after treating n as a symbolic parameter."
    )

    f3r = phi3(r)
    f6r = phi6(r)

    f3n = phi3(n)
    f6n = phi6(n)

    combinations = [
        ("Phi3(r),Phi3(n)", f3r, f3n),
        ("Phi3(r),Phi6(n)", f3r, f6n),
        ("Phi6(r),Phi3(n)", f6r, f3n),
        ("Phi6(r),Phi6(n)", f6r, f6n),
    ]

    for name, a, b in combinations:

        print()
        print(
            name
        )

        res_r = sp.factor(
            sp.resultant(
                a,
                b,
                r
            )
        )

        res_n = sp.factor(
            sp.resultant(
                a,
                b,
                n
            )
        )

        print(
            "Resultant eliminating r:"
        )

        print(
            res_r
        )

        print(
            "Resultant eliminating n:"
        )

        print(
            res_n
        )


# =============================================================================
# SECTION 11
# SEARCH FOR SMALL-ROOT POLYNOMIALS
# =============================================================================

def search_small_root_relations():

    section(
        "11. SEARCH FOR SMALL-ROOT POLYNOMIAL RELATIONS"
    )

    print()
    print(
        "Search linear/quadratic transformations"
    )

    print(
        "r -> a*n + b"
    )

    print(
        "and"
    )

    print(
        "r -> a*u + b"
    )

    print(
        "that turn Phi3(r) or Phi6(r) into a polynomial"
    )

    print(
        "with an obvious factor involving n."
    )

    a, b = sp.symbols(
        "a b"
    )

    candidate_forms = [
        "a*n+b",
        "a*u+b",
    ]

    for form_name in candidate_forms:

        if form_name == "a*n+b":
            z = a*n + b
        else:
            z = a*u + b

        print()
        print(
            f"FORM = {form_name}"
        )

        for k, f in [
            (3, phi3),
            (6, phi6),
        ]:

            expr = sp.expand(
                f(z)
            )

            print()
            print(
                f"Phi_{k}({form_name}) ="
            )

            print(
                sp.factor(expr)
            )

    print()
    print(
        "This section is exploratory:"
    )

    print(
        "the goal is to discover transformations that"
    )

    print(
        "make the hidden-root condition algebraically visible."
    )


# =============================================================================
# SECTION 12
# MULTI-TARGET GCD PATTERN
# =============================================================================

def multitarget_gcd_patterns():

    section(
        "12. MULTI-TARGET GCD PATTERN SEARCH"
    )

    targets = [
        (2, 3),
        (5, 11),
        (17, 23),
        (31, 47),
        (41, 43),
        (53, 59),
        (61, 67),
        (71, 73),
        (79, 83),
        (89, 97),
        (101, 103),
        (127, 131),
        (137, 139),
        (149, 151),
        (157, 163),
    ]

    expressions = {
        "Phi3(n)": lambda n: phi3(n),
        "Phi6(n)": lambda n: phi6(n),
        "Phi3(n)*Phi6(n)":
            lambda n: phi3(n) * phi6(n),
        "n^2+n+1":
            lambda n: n*n + n + 1,
        "n^2-n+1":
            lambda n: n*n - n + 1,
        "n^3-1":
            lambda n: n**3 - 1,
        "n^3+1":
            lambda n: n**3 + 1,
    }

    print()
    print(
        "For each target and small r, compare gcd patterns."
    )

    for expr_name, expr_fn in expressions.items():

        print()
        print(
            f"EXPRESSION = {expr_name}"
        )

        all_gcds = []

        for p, q in targets:

            nv = p * q

            E = expr_fn(nv)

            for rv in range(
                1,
                11
            ):

                g = gcd(
                    E,
                    phi3(rv)
                )

                if g > 1:
                    all_gcds.append(
                        g
                    )

        if not all_gcds:

            print(
                "  no nontrivial gcds"
            )

        else:

            unique = sorted(
                set(all_gcds)
            )

            print(
                f"  nontrivial gcd values = {unique}"
            )

            print(
                f"  count = {len(all_gcds)}"
            )


# =============================================================================
# SECTION 13
# FINAL SUMMARY
# =============================================================================

def final_summary():

    section(
        "13. FINAL RESEARCH SUMMARY"
    )

    print()
    print(
        "The experiment investigated:"
    )

    print()
    print(
        "1. exact Phi3/Phi6 identities"
    )

    print(
        "2. symbolic resultants"
    )

    print(
        "3. gcd(Phi_i(r),Phi_j(n))"
    )

    print(
        "4. transformed arguments"
    )

    print(
        "5. order-theoretic interpretation"
    )

    print(
        "6. n-only gcd factor searches"
    )

    print(
        "7. small-root polynomial transformations"
    )

    print(
        "8. multi-target gcd patterns"
    )

    print()
    print(
        "The high-value outcome would be an identity of the form"
    )

    print()
    print(
        "    factor_of_n = gcd(n, G(n,r))"
    )

    print(
        "where G is constructed from cyclotomic structure"
    )

    print(
        "and r can be selected from n alone."
    )

    print()
    print(
        "Equally valuable would be an exact relation showing that"
    )

    print(
        "a gcd(Phi_i(r),Phi_j(n)) encodes a hidden divisor of n."
    )

    print()
    print(
        "NO CSV FILES WERE CREATED."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print()
    print("=" * 78)
    print(
        "CYCLOTOMIC RESULTANT / GCD / ORDER-STRUCTURE EXPERIMENT"
    )
    print("=" * 78)

    test_exact_identities()

    test_resultants()

    test_symbolic_gcds()

    test_numeric_gcd_table()

    factor_gcd_patterns()

    substitution_search()

    numerical_substitution_gcds()

    order_interpretation()

    resultant_candidate_search()

    dual_n_resultants()

    search_small_root_relations()

    multitarget_gcd_patterns()

    final_summary()

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()

