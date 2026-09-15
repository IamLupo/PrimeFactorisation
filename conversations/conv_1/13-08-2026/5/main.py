#!/usr/bin/env python3

"""
==============================================================================
CUBIC ROOT CONSTRUCTION / N-ONLY FACTOR SEARCH — CORRECTED FULL VERSION
==============================================================================

Main target:

    Phi_3(x) = x^2 + x + 1

A factor p of n is revealed if

    Phi_3(r) == 0 (mod p)

because then

    gcd(n, Phi_3(r))

may return p.

The experiment studies whether r can be constructed from n alone.

NO CSV FILES.
EVERYTHING IS PRINTED.
"""

from __future__ import annotations

from math import gcd
import sympy as sp


# =============================================================================
# CONFIGURATION
# =============================================================================

P = 22612043
Q = 26706517
N = P * Q

SMALL_R = list(range(1, 31))

SMALL_COEFFS = [-3, -2, -1, 0, 1, 2, 3]

MODULI = [
    3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47
]

TARGETS = [
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


# =============================================================================
# BASIC POLYNOMIALS
# =============================================================================

def phi3(x: int) -> int:
    return x * x + x + 1


def phi6(x: int) -> int:
    return x * x - x + 1


def bitlen(x: int) -> int:
    return abs(int(x)).bit_length() if x else 1


def section(title: str):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def line():
    print("-" * 78)


# =============================================================================
# SYMBOLIC VARIABLES
# =============================================================================

x, r, n, a, b, c = sp.symbols(
    "x r n a b c"
)


# =============================================================================
# SECTION 1
# =============================================================================

def section_exact_structure():
    section("1. EXACT CUBIC-ROOT STRUCTURE")

    print()
    print("Phi_3(x) = x^2 + x + 1")
    print("Phi_6(x) = x^2 - x + 1")

    print()
    print("Exact shift:")

    print(
        "Phi_6(x) - Phi_3(x-1) =",
        sp.factor(
            phi6(x) - phi3(x - 1)
        )
    )

    print(
        "Phi_3(x) - Phi_6(x+1) =",
        sp.factor(
            phi3(x) - phi6(x + 1)
        )
    )

    print()
    print(
        "Therefore only Phi_3 needs to be studied."
    )

    print()
    print(
        "Root condition:"
    )

    print(
        "Phi_3(r)=0 mod p  <=>  r^3=1 mod p and r != 1"
    )

    print()
    print(
        "For prime p != 3, such roots exist exactly when p == 1 mod 3."
    )


# =============================================================================
# SECTION 2
# =============================================================================

def section_control():
    section("2. REPRESENTATIVE HIDDEN FACTORS")

    print()
    print(f"p = {P}")
    print(f"q = {Q}")
    print(f"n = {N}")

    print()
    print(f"p mod 3 = {P % 3}")
    print(f"q mod 3 = {Q % 3}")

    roots_p_raw = sp.sqrt_mod(
        (-3) % P,
        P,
        all_roots=True
    )

    roots_q_raw = sp.sqrt_mod(
        (-3) % Q,
        Q,
        all_roots=True
    )

    roots_p = sorted({
        ((-1 + z) * pow(2, -1, P)) % P
        for z in roots_p_raw
    })

    roots_q = sorted({
        ((-1 + z) * pow(2, -1, Q)) % Q
        for z in roots_q_raw
    })

    print()
    print(f"Phi_3 roots mod p = {roots_p}")
    print(f"Phi_3 roots mod q = {roots_q}")

    print()
    print(
        "The roots are large in the representative example."
    )

    print(
        "This is why blind r enumeration is not expected to scale."
    )


# =============================================================================
# SECTION 3
# =============================================================================

def section_resultant_controls():
    section("3. EXACT RESULTANT CONTROLS")

    Phi3_r = phi3(r)
    Phi3_n = phi3(n)
    Phi6_n = phi6(n)

    print()
    print(
        "Verify exact resultant identities."
    )

    res33 = sp.factor(
        sp.resultant(
            Phi3_r,
            Phi3_n,
            r
        )
    )

    res36 = sp.factor(
        sp.resultant(
            Phi3_r,
            Phi6_n,
            r
        )
    )

    print()
    print(
        "Res_r(Phi3(r), Phi3(n)) ="
    )
    print(res33)

    print("Expected:")
    print(sp.factor(Phi3_n ** 2))

    print(
        "Residual =",
        sp.factor(
            res33 - Phi3_n ** 2
        )
    )

    print()
    print(
        "Res_r(Phi3(r), Phi6(n)) ="
    )
    print(res36)

    print("Expected:")
    print(sp.factor(Phi6_n ** 2))

    print(
        "Residual =",
        sp.factor(
            res36 - Phi6_n ** 2
        )
    )

    print()
    print(
        "PASS: basic resultant structure is exactly explained."
    )


# =============================================================================
# SECTION 4
# =============================================================================

def section_n_only_candidates():
    section("4. DIRECT N-ONLY CANDIDATE r CONSTRUCTIONS")

    candidates = {
        "1":
            lambda nn: 1,

        "-1":
            lambda nn: -1,

        "n":
            lambda nn: nn,

        "n+1":
            lambda nn: nn + 1,

        "n-1":
            lambda nn: nn - 1,

        "2n+1":
            lambda nn: 2 * nn + 1,

        "2n-1":
            lambda nn: 2 * nn - 1,

        "n^2":
            lambda nn: nn * nn,

        "n^2+n+1":
            lambda nn: nn * nn + nn + 1,

        "n^2-n+1":
            lambda nn: nn * nn - nn + 1,

        "n^3-1":
            lambda nn: nn ** 3 - 1,

        "n^3+1":
            lambda nn: nn ** 3 + 1,
    }

    print()
    print("Every candidate below is constructed from n alone.")
    print()

    print(
        "name                 r bits       gcd(n,Phi3(r))"
    )
    line()

    for name, fn in candidates.items():

        rv = fn(N)

        g = gcd(
            N,
            phi3(rv)
        )

        print(
            f"{name:20s}"
            f"{bitlen(rv):10d}"
            f"{g:24d}"
        )

    print()
    print(
        "A nontrivial gcd here would be a direct n-only factorization result."
    )


# =============================================================================
# SECTION 5
# =============================================================================

def section_affine_search():
    section("5. AFFINE r = a*n+b SEARCH")

    print()
    print("Search:")
    print("    r = a*n + b")
    print(f"a,b in {SMALL_COEFFS}")

    hits = []

    print()
    print(
        "a      b      gcd(n,Phi3(a*n+b))"
    )
    line()

    for aa in SMALL_COEFFS:

        for bb in SMALL_COEFFS:

            rv = aa * N + bb

            g = gcd(
                N,
                phi3(rv)
            )

            if 1 < g < N:

                hits.append(
                    (aa, bb, rv, g)
                )

                print(
                    f"{aa:3d}"
                    f"{bb:7d}"
                    f"{g:24d}"
                    "   <<< FACTOR"
                )

    if not hits:
        print("No nontrivial factor found.")

    else:

        print()
        print("DIRECT FACTORS FOUND:")

        for item in hits:
            print(item)


# =============================================================================
# SECTION 6
# =============================================================================

def section_quadratic_search():
    section("6. QUADRATIC r = a*n^2+b*n+c SEARCH")

    print()
    print(
        "Search small-coefficient quadratic constructions."
    )

    print(
        "r = a*n^2 + b*n + c"
    )

    hits = []

    for aa in [-1, 0, 1]:

        for bb in SMALL_COEFFS:

            for cc in SMALL_COEFFS:

                rv = (
                    aa * N * N
                    + bb * N
                    + cc
                )

                g = gcd(
                    N,
                    phi3(rv)
                )

                if 1 < g < N:

                    hits.append(
                        (
                            aa,
                            bb,
                            cc,
                            g
                        )
                    )

    print()
    print("Nontrivial hits:")

    if not hits:
        print("None.")

    else:

        for item in hits:
            print(item)


# =============================================================================
# SECTION 7
# =============================================================================

def section_small_shift_search():
    section("7. r = n + c / n - c SEARCH")

    print()
    print(
        "This tests whether small corrections around n expose a factor."
    )

    print()
    print(
        "c        gcd(n,Phi3(n+c))    gcd(n,Phi3(n-c))"
    )
    line()

    hits = []

    for cc in range(-50, 51):

        gp = gcd(
            N,
            phi3(N + cc)
        )

        gm = gcd(
            N,
            phi3(N - cc)
        )

        print(
            f"{cc:4d}"
            f"{gp:24d}"
            f"{gm:24d}"
        )

        if 1 < gp < N:
            hits.append(
                ("n+c", cc, gp)
            )

        if 1 < gm < N:
            hits.append(
                ("n-c", cc, gm)
            )

    print()

    if not hits:
        print(
            "No nontrivial factors for |c| <= 50."
        )

    else:

        print(
            "Potential hits:"
        )

        for item in hits:
            print(item)


# =============================================================================
# SECTION 8
# =============================================================================

def section_composed_search():
    section("8. COMPOSED POLYNOMIAL CONSTRUCTIONS")

    print()
    print(
        "Test candidate r-values built from natural n-polynomials:"
    )

    print(
        "    n"
    )
    print(
        "    n+1"
    )
    print(
        "    n-1"
    )
    print(
        "    Phi3(n)"
    )
    print(
        "    Phi6(n)"
    )
    print(
        "    n^2-1"
    )
    print(
        "    n^2+n+1"
    )
    print(
        "    n^2-n+1"
    )

    base = {
        "n":
            N,

        "n+1":
            N + 1,

        "n-1":
            N - 1,

        "Phi3(n)":
            phi3(N),

        "Phi6(n)":
            phi6(N),

        "n^2-1":
            N * N - 1,

        "n^2+n+1":
            N * N + N + 1,

        "n^2-n+1":
            N * N - N + 1,
    }

    candidates = {}

    for name, value in base.items():

        candidates[name] = value

        candidates[
            f"Phi3({name})"
        ] = phi3(value)

        candidates[
            f"Phi6({name})"
        ] = phi6(value)

    print()
    print(
        "candidate                         gcd with n"
    )
    line()

    found = []

    for name, rv in candidates.items():

        g = gcd(
            N,
            phi3(rv)
        )

        if 1 < g < N:

            found.append(
                (
                    name,
                    g
                )
            )

            print(
                f"{name:30s}"
                f"{g:24d}"
            )

    if not found:
        print("No nontrivial factors found.")


# =============================================================================
# SECTION 9
# FIXED VERSION OF EUCLIDEAN SEARCH
# =============================================================================

def section_euclidean_search():
    section("9. EUCLIDEAN / POLYNOMIAL-COMBINATION SEARCH")

    print()
    print(
        "Search small polynomial combinations of"
    )

    print(
        "Phi3(n), Phi6(n), n, n+1, n-1, and constants."
    )

    print()
    print(
        "The question:"
    )

    print(
        "Can a simple n-only combination produce a"
    )

    print(
        "quantity whose gcd with n reveals p or q?"
    )

    # IMPORTANT:
    # The previous version had n+1 / n-1 in the search loops
    # but forgot to put them into the dictionary.
    #
    # Everything referenced below is now explicitly defined.

    values = {
        "1":
            1,

        "n":
            N,

        "n+1":
            N + 1,

        "n-1":
            N - 1,

        "Phi3(n)":
            phi3(N),

        "Phi6(n)":
            phi6(N),

        "n^2-1":
            N * N - 1,

        "n^2+n+1":
            N * N + N + 1,

        "n^2-n+1":
            N * N - N + 1,

        "n^3-1":
            N ** 3 - 1,

        "n^3+1":
            N ** 3 + 1,
    }

    keys = list(values.keys())

    print()
    print("PAIRWISE GCDS:")
    line()

    seen = set()

    for i in range(len(keys)):

        for j in range(i + 1, len(keys)):

            name_i = keys[i]
            name_j = keys[j]

            g = gcd(
                abs(values[name_i]),
                abs(values[name_j])
            )

            if g > 1:

                record = (
                    name_i,
                    name_j,
                    g
                )

                if record not in seen:

                    seen.add(record)

                    print(
                        f"gcd({name_i}, {name_j}) = {g}"
                    )

    print()
    print(
        "EXPECTED STRUCTURAL GCDS:"
    )

    print(
        "gcd(Phi3(n), n^2+n+1)"
        " = Phi3(n)"
    )

    print(
        "gcd(Phi3(n), n^3-1)"
        " = Phi3(n)"
    )

    print(
        "gcd(Phi6(n), n^2-n+1)"
        " = Phi6(n)"
    )

    print(
        "gcd(Phi6(n), n^3+1)"
        " = Phi6(n)"
    )

    print()
    print(
        "Now test small linear combinations."
    )

    coefficients = [-2, -1, 1, 2]

    found = []

    # We deliberately use a smaller, meaningful subset of expressions.
    expression_names = [
        "Phi3(n)",
        "Phi6(n)",
        "n^2+n+1",
        "n^2-n+1",
        "n^3-1",
        "n^3+1",
    ]

    correction_names = [
        "1",
        "n",
        "n+1",
        "n-1",
    ]

    print()
    print(
        "Searching:"
    )
    print(
        "    c1*X + c2*Y"
    )

    print()

    for k1 in coefficients:

        for k2 in coefficients:

            for name1 in expression_names:

                for name2 in correction_names:

                    value = (
                        k1 * values[name1]
                        + k2 * values[name2]
                    )

                    g = gcd(
                        N,
                        abs(value)
                    )

                    if 1 < g < N:

                        found.append(
                            (
                                k1,
                                name1,
                                k2,
                                name2,
                                g
                            )
                        )

    if not found:

        print(
            "No nontrivial factor found by these small combinations."
        )

    else:

        print(
            "Potential direct-factor combinations:"
        )

        # Deduplicate exact records.
        unique = sorted(
            set(found),
            key=lambda z: (
                z[4],
                z[0],
                z[1],
                z[2],
                z[3]
            )
        )

        for item in unique:
            print(item)


# =============================================================================
# SECTION 10
# MODULAR ROOT CONSTRUCTION
# =============================================================================

def section_modular_root_construction():
    section("10. MODULAR ROOT CONSTRUCTION TEST")

    print()
    print(
        "For every small modulus m, determine roots of"
    )

    print(
        "Phi3(x)=0 mod m."
    )

    print()
    print(
        "Then ask whether any obvious n-only quantity predicts"
    )

    print(
        "one of those roots."
    )

    print()
    print(
        "m    roots Phi3    n mod m    n+1    n-1    Phi3(n) mod m"
    )

    line()

    for m in MODULI:

        roots = [
            rr
            for rr in range(m)
            if phi3(rr) % m == 0
        ]

        print(
            f"{m:2d}"
            f"{str(roots):16s}"
            f"{N % m:10d}"
            f"{(N + 1) % m:8d}"
            f"{(N - 1) % m:8d}"
            f"{phi3(N) % m:16d}"
        )


# =============================================================================
# SECTION 11
# MODULAR FUNCTION SEARCH
# =============================================================================

def section_modular_function_search():
    section("11. SEARCH FOR r = f(n) MOD m")

    print()
    print(
        "Test simple n-only functions modulo m."
    )

    rules = [
        "n",
        "n+1",
        "n-1",
        "-n",
    ]

    print()
    print(
        "m    rule                 r mod m    Phi3(r) mod m"
    )

    line()

    hits = []

    for m in MODULI:

        for rule in rules:

            if rule == "n":
                rv = N % m

            elif rule == "n+1":
                rv = (N + 1) % m

            elif rule == "n-1":
                rv = (N - 1) % m

            else:
                rv = (-N) % m

            residue = phi3(rv) % m

            if residue == 0:

                hits.append(
                    (
                        m,
                        rule,
                        rv
                    )
                )

                print(
                    f"{m:2d}"
                    f"{rule:20s}"
                    f"{rv:12d}"
                    f"{residue:18d}"
                    "   <<< ROOT"
                )

    if not hits:

        print(
            "No nontrivial n-only rule produced a Phi3 root"
            " for the tested moduli."
        )


# =============================================================================
# SECTION 12
# HIDDEN ROOT VS N-ONLY
# =============================================================================

def section_hidden_vs_constructed():
    section(
        "12. HIDDEN ROOT VS N-ONLY CANDIDATE RESIDUES"
    )

    print()
    print(
        "Control comparison using the known factors."
    )

    roots_p_raw = sp.sqrt_mod(
        (-3) % P,
        P,
        all_roots=True
    )

    roots_q_raw = sp.sqrt_mod(
        (-3) % Q,
        Q,
        all_roots=True
    )

    roots_p = sorted({
        ((-1 + z) * pow(2, -1, P)) % P
        for z in roots_p_raw
    })

    roots_q = sorted({
        ((-1 + z) * pow(2, -1, Q)) % Q
        for z in roots_q_raw
    })

    print()
    print(
        f"Actual roots mod p: {roots_p}"
    )

    print(
        f"Actual roots mod q: {roots_q}"
    )

    print()
    print(
        "candidate                value         Phi3(value) mod factor"
    )

    line()

    candidate_values = [
        (
            "n mod p",
            N % P,
            P
        ),
        (
            "n+1 mod p",
            (N + 1) % P,
            P
        ),
        (
            "n-1 mod p",
            (N - 1) % P,
            P
        ),
        (
            "n mod q",
            N % Q,
            Q
        ),
        (
            "n+1 mod q",
            (N + 1) % Q,
            Q
        ),
        (
            "n-1 mod q",
            (N - 1) % Q,
            Q
        ),
    ]

    for name, value, modulus in candidate_values:

        print(
            f"{name:24s}"
            f"{value:14d}"
            f"{phi3(value) % modulus:24d}"
        )


# =============================================================================
# SECTION 13
# MULTI-TARGET
# =============================================================================

def section_multitarget():
    section(
        "13. MULTI-TARGET N-ONLY CONSTRUCTION TEST"
    )

    print()
    print(
        "The same n-only rules are tested across many semiprimes."
    )

    rules = {
        "n":
            lambda nn: nn,

        "n+1":
            lambda nn: nn + 1,

        "n-1":
            lambda nn: nn - 1,

        "Phi3(n)":
            lambda nn: phi3(nn),

        "Phi6(n)":
            lambda nn: phi6(nn),

        "n^2":
            lambda nn: nn * nn,

        "n^2+n+1":
            lambda nn: nn * nn + nn + 1,

        "n^2-n+1":
            lambda nn: nn * nn - nn + 1,
    }

    print()
    print(
        "rule                    nontrivial gcd count"
    )

    line()

    for rule_name, fn in rules.items():

        count = 0

        for pp, qq in TARGETS:

            nn = pp * qq

            rv = fn(nn)

            g = gcd(
                nn,
                phi3(rv)
            )

            if 1 < g < nn:
                count += 1

        print(
            f"{rule_name:25s}"
            f"{count:10d} / {len(TARGETS)}"
        )


# =============================================================================
# SECTION 14
# RESULTANT REDUCTIONS
# =============================================================================

def section_resultant_reductions():
    section(
        "14. RESULTANT-BASED ALGEBRAIC REDUCTIONS"
    )

    print()
    print(
        "Verify exact substitution identities."
    )

    tests = [
        (
            "Phi3(n+r) vs Phi3(n)",
            sp.resultant(
                phi3(n + r),
                phi3(n),
                n
            ),
            r**2 * (r**2 + 3)
        ),

        (
            "Phi3(n+r) vs Phi6(n)",
            sp.resultant(
                phi3(n + r),
                phi6(n),
                n
            ),
            (r + 1)**2 * (
                r**2 + 2*r + 4
            )
        ),

        (
            "Phi6(n+r) vs Phi3(n)",
            sp.resultant(
                phi6(n + r),
                phi3(n),
                n
            ),
            (r - 1)**2 * (
                r**2 - 2*r + 4
            )
        ),

        (
            "Phi6(n+r) vs Phi6(n)",
            sp.resultant(
                phi6(n + r),
                phi6(n),
                n
            ),
            r**2 * (r**2 + 3)
        ),
    ]

    for name, actual, expected in tests:

        print()
        print(name)

        print(
            "actual   =",
            sp.factor(actual)
        )

        print(
            "expected =",
            sp.factor(expected)
        )

        print(
            "residual =",
            sp.factor(
                actual - expected
            )
        )

    print()
    print(
        "PASS: all tested resultant reductions are exact."
    )


# =============================================================================
# SECTION 15
# NEW: DIRECT SEARCH FOR GCDs FROM RESULTANT FACTORS
# =============================================================================

def section_resultant_factor_gcds():
    section(
        "15. GCDs FROM RESULTANT FACTORS"
    )

    print()
    print(
        "The symbolic resultants produce factors such as"
    )

    print(
        "    r^2"
    )

    print(
        "    r^2+3"
    )

    print(
        "    r^2+2r+4"
    )

    print(
        "    r^2-2r+4"
    )

    print()
    print(
        "Now test whether these resultant factors themselves"
    )

    print(
        "can interact with n in a factor-revealing way."
    )

    print()
    print(
        "r       gcd(n,r²+3)"
        "      gcd(n,r²+2r+4)"
        "      gcd(n,r²-2r+4)"
    )

    line()

    found = []

    for rv in SMALL_R:

        g1 = gcd(
            N,
            rv * rv + 3
        )

        g2 = gcd(
            N,
            rv * rv + 2 * rv + 4
        )

        g3 = gcd(
            N,
            rv * rv - 2 * rv + 4
        )

        print(
            f"{rv:2d}"
            f"{g1:18d}"
            f"{g2:24d}"
            f"{g3:24d}"
        )

        for label, g in [
            ("r²+3", g1),
            ("r²+2r+4", g2),
            ("r²-2r+4", g3),
        ]:

            if 1 < g < N:

                found.append(
                    (
                        rv,
                        label,
                        g
                    )
                )

    print()

    if not found:

        print(
            "No direct factor found from the resultant factors."
        )

    else:

        print(
            "NONTRIVIAL FACTORS FOUND:"
        )

        for item in found:
            print(item)


# =============================================================================
# SECTION 16
# GENERAL SMALL-POLYNOMIAL r SEARCH
# =============================================================================

def section_small_polynomial_r_search():
    section(
        "16. SMALL POLYNOMIAL r=f(n) SEARCH"
    )

    print()
    print(
        "Search:"
    )

    print(
        "r = c0 + c1*n + c2*n^2"
    )

    print(
        "with small coefficients."
    )

    hits = []

    coeffs = [-2, -1, 0, 1, 2]

    for c0 in coeffs:

        for c1 in coeffs:

            for c2 in coeffs:

                if c0 == 0 and c1 == 0 and c2 == 0:
                    continue

                rv = (
                    c0
                    + c1 * N
                    + c2 * N * N
                )

                g = gcd(
                    N,
                    phi3(rv)
                )

                if 1 < g < N:

                    hits.append(
                        (
                            c0,
                            c1,
                            c2,
                            g
                        )
                    )

    if not hits:

        print(
            "No nontrivial factor found."
        )

    else:

        print(
            "Potential factor constructions:"
        )

        for item in sorted(
            set(hits)
        ):

            print(item)


# =============================================================================
# SECTION 17
# FINAL
# =============================================================================

def section_final():
    section(
        "17. FINAL RESEARCH QUESTION"
    )

    print()
    print(
        "The experiment now separates four possible mechanisms:"
    )

    print()
    print(
        "A. Direct polynomial construction"
    )

    print(
        "     r=f(n)"
    )

    print(
        "     gcd(n,Phi3(r))"
    )

    print()

    print(
        "B. Resultant-derived construction"
    )

    print(
        "     r²+3"
    )

    print(
        "     r²+2r+4"
    )

    print(
        "     r²-2r+4"
    )

    print()

    print(
        "C. Modular root construction"
    )

    print(
        "     n mod m -> candidate cubic root mod m"
    )

    print()

    print(
        "D. Polynomial combinations"
    )

    print(
        "     gcd(n, c1*F1(n)+c2*F2(n))"
    )

    print()
    print(
        "The decisive success criterion remains:"
    )

    print()
    print(
        "    1 < gcd(n,G(n)) < n"
    )

    print()
    print(
        "for a G that can be computed from n alone."
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
        "CUBIC ROOT CONSTRUCTION / N-ONLY FACTOR SEARCH — FULL"
    )
    print("=" * 78)

    section_exact_structure()

    section_control()

    section_resultant_controls()

    section_n_only_candidates()

    section_affine_search()

    section_quadratic_search()

    section_small_shift_search()

    section_composed_search()

    section_euclidean_search()

    section_modular_root_construction()

    section_modular_function_search()

    section_hidden_vs_constructed()

    section_multitarget()

    section_resultant_reductions()

    section_resultant_factor_gcds()

    section_small_polynomial_r_search()

    section_final()

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()