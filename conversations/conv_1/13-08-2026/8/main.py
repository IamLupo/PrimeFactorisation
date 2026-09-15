#!/usr/bin/env python3

"""
==============================================================================
CUBIC ROOT GENERATION BY MODULAR EXPONENTIATION
==============================================================================

NEW DIRECTION
-------------

The previous experiments establish a strong obstruction:

For any polynomial G,

    G(n) == G(0) (mod n)

therefore

    gcd(n, G(n)) = gcd(n, G(0)).

So pure polynomial manipulation of n cannot reveal a hidden factor.

We therefore move to a genuinely different operation:

    modular exponentiation.

The cubic-root target is

    x^2 + x + 1 == 0 (mod p)

which is equivalent to

    x^3 == 1 (mod p), x != 1.

Instead of guessing x, we generate

    x = a^E mod n

for many bases a and exponents E.

If x has order 3 modulo one hidden factor but not the other, then

    gcd(n, x^2+x+1)

may reveal that factor.

This experiment tests:

1. powers a^E for many E
2. exponents related to 3
3. exponents related to lcm(1,...,B)
4. smooth exponents
5. repeated powering
6. whether the gcd success is correlated with p-1 or q-1
7. whether cubic-root generation is easier than directly finding
   a cubic root
8. comparison with ordinary Pollard p-1 style behavior

NO CSV FILES.
EVERYTHING IS PRINTED.
"""


from math import gcd, lcm
import sympy as sp


# =============================================================================
# REPRESENTATIVE TARGET
# =============================================================================

P = 22612043
Q = 26706517
N = P * Q


# =============================================================================
# CONFIGURATION
# =============================================================================

BASES = list(range(2, 31))

# Direct exponent families.
EXPONENTS = [
    1,
    2,
    3,
    4,
    5,
    6,
    9,
    12,
    15,
    18,
    21,
    24,
    27,
    30,
    36,
    45,
    60,
    90,
    120,
    180,
    240,
    360,
    720,
]

SMOOTH_BOUNDS = [
    5,
    7,
    10,
    12,
    15,
    20,
    25,
    30,
    40,
]


# =============================================================================
# HELPERS
# =============================================================================

def phi3(x):
    return x * x + x + 1


def bitlen(x):
    return abs(int(x)).bit_length() if x else 1


def section(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def line():
    print("-" * 78)


def nontrivial_factor(g):
    if 1 < g < N:
        return g
    return None


# =============================================================================
# SECTION 1
# POLYNOMIAL OBSTRUCTION
# =============================================================================

def test_polynomial_obstruction():

    section(
        "1. POLYNOMIAL N-ONLY OBSTRUCTION"
    )

    print()
    print(
        "For every integer polynomial G:"
    )

    print()
    print(
        "    G(n) == G(0) (mod n)"
    )

    print()
    print(
        "Therefore:"
    )

    print()
    print(
        "    gcd(n,G(n)) = gcd(n,G(0))"
    )

    print()
    print(
        "Examples:"
    )

    examples = {
        "n+1": N + 1,
        "n-1": N - 1,
        "n^2+n+1": N*N + N + 1,
        "n^2-n+1": N*N - N + 1,
        "n^3-1": N**3 - 1,
        "n^3+1": N**3 + 1,
    }

    for name, value in examples.items():

        print(
            f"{name:15s}"
            f" value mod n = {value % N}"
            f" gcd = {gcd(N,value)}"
        )

    print()
    print(
        "CONCLUSION:"
    )

    print(
        "Pure polynomial substitution is mathematically exhausted."
    )


# =============================================================================
# SECTION 2
# FACTOR ORDER DATA
# =============================================================================

def factor_order_data():

    section(
        "2. HIDDEN FACTOR ORDER DATA"
    )

    print()
    print(
        f"p = {P}"
    )

    print(
        f"q = {Q}"
    )

    print(
        f"p-1 = {P-1}"
    )

    print(
        f"q-1 = {Q-1}"
    )

    print()

    print(
        "Factorizations:"
    )

    print(
        f"p-1 = {sp.factorint(P-1)}"
    )

    print(
        f"q-1 = {sp.factorint(Q-1)}"
    )

    print()

    print(
        "3-adic valuations:"
    )

    def valuation_3(x):
        count = 0
        while x % 3 == 0:
            count += 1
            x //= 3
        return count

    print(
        f"v3(p-1) = {valuation_3(P-1)}"
    )

    print(
        f"v3(q-1) = {valuation_3(Q-1)}"
    )

    print()
    print(
        "This tells us whether order-3 elements are structurally"
    )

    print(
        "available in each hidden factor."
    )


# =============================================================================
# SECTION 3
# DIRECT POWER SEARCH
# =============================================================================

def direct_power_search():

    section(
        "3. DIRECT MODULAR POWER SEARCH"
    )

    print()
    print(
        "For each base a and exponent E:"
    )

    print()
    print(
        "    x = a^E mod n"
    )

    print(
        "    g = gcd(n, x^2+x+1)"
    )

    print()
    print(
        "A nontrivial g is a successful cubic-root factor discovery."
    )

    print()
    print(
        "a     E       x bits       gcd"
    )

    line()

    hits = []

    for a in BASES:

        if gcd(a,N) != 1:
            continue

        for E in EXPONENTS:

            x = pow(a,E,N)

            g = gcd(
                N,
                phi3(x)
            )

            print(
                f"{a:2d}"
                f"{E:8d}"
                f"{bitlen(x):12d}"
                f"{g:24d}"
            )

            if 1 < g < N:

                hits.append(
                    (
                        a,
                        E,
                        x,
                        g
                    )
                )

    print()

    if not hits:

        print(
            "No nontrivial factor found."
        )

    else:

        print(
            "SUCCESSFUL CUBIC-ROOT FACTOR HITS:"
        )

        for a, E, x, g in hits:

            print()
            print(
                f"base={a}, exponent={E}"
            )

            print(
                f"x={x}"
            )

            print(
                f"gcd={g}"
            )


# =============================================================================
# SECTION 4
# TEST MODULO HIDDEN FACTORS
# =============================================================================

def order_control():

    section(
        "4. CONTROL: WHAT HAPPENS MODULO p AND q?"
    )

    print()
    print(
        "For selected powers, inspect:"
    )

    print(
        "    x_p = a^E mod p"
    )

    print(
        "    x_q = a^E mod q"
    )

    print(
        "and test Phi3(x_p), Phi3(x_q)."
    )

    tests = [
        (2, 3),
        (2, 6),
        (2, 9),
        (2, 12),
        (2, 18),
        (2, 21),
        (2, 30),
        (3, 3),
        (5, 3),
        (7, 3),
    ]

    print()
    print(
        "a    E       x mod p   Phi3(x) mod p"
        "      x mod q   Phi3(x) mod q"
    )

    line()

    for a, E in tests:

        xp = pow(a,E,P)
        xq = pow(a,E,Q)

        print(
            f"{a:2d}"
            f"{E:5d}"
            f"{xp:12d}"
            f"{phi3(xp) % P:18d}"
            f"{xq:12d}"
            f"{phi3(xq) % Q:18d}"
        )


# =============================================================================
# SECTION 5
# EXPONENTS WITH A FACTOR OF 3
# =============================================================================

def exponent_3_structure():

    section(
        "5. EXPONENTS CONTAINING POWERS OF 3"
    )

    print()
    print(
        "Search:"
    )

    print(
        "    E = 3^k"
    )

    print(
        "and"
    )

    print(
        "    E = 3^k * t"
    )

    print(
        "The aim is to see whether repeated cubing naturally"
    )

    print(
        "creates order-3 elements modulo one hidden factor."
    )

    hits = []

    for a in BASES:

        if gcd(a,N) != 1:
            continue

        for k in range(0, 15):

            E = 3**k

            x = pow(a,E,N)

            g = gcd(
                N,
                phi3(x)
            )

            if 1 < g < N:

                hits.append(
                    (
                        a,
                        E,
                        g
                    )
                )

            print(
                f"a={a:2d}"
                f" k={k:2d}"
                f" E={E:8d}"
                f" gcd={g}"
            )

    print()

    if not hits:

        print(
            "No cubic-root factor hit from pure powers of 3."
        )

    else:

        print(
            "HITS:"
        )

        for hit in hits:
            print(hit)


# =============================================================================
# SECTION 6
# SMOOTH EXPONENT SEARCH
# =============================================================================

def smooth_exponent(B):

    E = 1

    for prime in sp.primerange(2, B + 1):

        power = prime

        while power * prime <= B:

            power *= prime

        E *= power

    return int(E)


def smooth_search():

    section(
        "6. SMOOTH-EXPONENT SEARCH"
    )

    print()
    print(
        "For each B:"
    )

    print(
        "    E = lcm(1,2,...,B)"
    )

    print(
        "Then test a^E modulo n."
    )

    print()
    print(
        "B       E bits      base      gcd"
    )

    line()

    hits = []

    for B in SMOOTH_BOUNDS:

        E = smooth_exponent(B)

        for a in BASES:

            if gcd(a,N) != 1:
                continue

            x = pow(a,E,N)

            g = gcd(
                N,
                phi3(x)
            )

            print(
                f"{B:2d}"
                f"{bitlen(E):12d}"
                f"{a:10d}"
                f"{g:24d}"
            )

            if 1 < g < N:

                hits.append(
                    (
                        B,
                        E,
                        a,
                        x,
                        g
                    )
                )

    print()

    if not hits:

        print(
            "No factor found using the smooth exponent search."
        )

    else:

        print(
            "SUCCESSFUL HITS:"
        )

        for item in hits:
            print(item)


# =============================================================================
# SECTION 7
# POLLARD-STYLE COMPARISON
# =============================================================================

def pollard_p1_comparison():

    section(
        "7. POLLARD-p-1 STYLE CONTROL"
    )

    print()
    print(
        "For the same smooth exponents, compare:"
    )

    print(
        "    gcd(n, a^E - 1)"
    )

    print(
        "against:"
    )

    print(
        "    gcd(n, (a^E)^2 + a^E + 1)"
    )

    print()
    print(
        "This tells us whether the cubic-root target gives"
    )

    print(
        "anything beyond ordinary order-finding."
    )

    print()
    print(
        "B    a     gcd(a^E-1)       gcd(Phi3(a^E))"
    )

    line()

    for B in SMOOTH_BOUNDS:

        E = smooth_exponent(B)

        for a in BASES:

            if gcd(a,N) != 1:
                continue

            x = pow(a,E,N)

            g1 = gcd(
                N,
                x - 1
            )

            g2 = gcd(
                N,
                phi3(x)
            )

            if (
                1 < g1 < N
                or
                1 < g2 < N
            ):

                print(
                    f"{B:2d}"
                    f"{a:5d}"
                    f"{g1:20d}"
                    f"{g2:24d}"
                )


# =============================================================================
# SECTION 8
# ROOT GENERATION VIA GROUP ORDER
# =============================================================================

def group_order_experiment():

    section(
        "8. TRY TO FORCE AN ELEMENT OF ORDER 3"
    )

    print()
    print(
        "If we knew p-1 = 3*m, then"
    )

    print(
        "    a^m"
    )

    print(
        "would have order dividing 3 modulo p."
    )

    print()
    print(
        "We do not know p-1."
    )

    print(
        "So we test candidate exponents generated from"
    )

    print(
        "small smooth numbers."
    )

    print()
    print(
        "m            base    gcd(Phi3(a^m), n)"
    )

    line()

    hits = []

    candidate_m = []

    for B in range(2, 31):

        candidate_m.append(
            smooth_exponent(B)
        )

    candidate_m = sorted(
        set(candidate_m)
    )

    for m in candidate_m:

        for a in BASES:

            if gcd(a,N) != 1:
                continue

            x = pow(
                a,
                m,
                N
            )

            g = gcd(
                N,
                phi3(x)
            )

            if 1 < g < N:

                hits.append(
                    (
                        m,
                        a,
                        g
                    )
                )

                print(
                    f"{m:18d}"
                    f"{a:8d}"
                    f"{g:24d}"
                )

    if not hits:

        print(
            "No forced order-3 factor was found."
        )


# =============================================================================
# SECTION 9
# ROOT OF UNITY CHECK
# =============================================================================

def root_of_unity_check():

    section(
        "9. ROOT-OF-UNITY CHECK"
    )

    print()
    print(
        "Whenever x is found, check:"
    )

    print(
        "    x^3 mod p"
    )

    print(
        "    x^3 mod q"
    )

    print(
        "and compare with the cubic-root condition."
    )

    # Use some of the direct powers.
    tests = [
        (2, 3),
        (2, 6),
        (2, 9),
        (3, 6),
        (5, 6),
        (7, 12),
        (11, 12),
    ]

    print()
    print(
        "a    E    x_p^3 mod p    x_q^3 mod q"
    )

    line()

    for a, E in tests:

        xp = pow(a,E,P)
        xq = pow(a,E,Q)

        print(
            f"{a:2d}"
            f"{E:5d}"
            f"{pow(xp,3,P):18d}"
            f"{pow(xq,3,Q):18d}"
        )


# =============================================================================
# SECTION 10
# MULTI-TARGET SCALING
# =============================================================================

def multi_target_scaling():

    section(
        "10. MULTI-TARGET SCALING"
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
        (179, 181),
        (191, 193),
    ]

    B = 20
    E = smooth_exponent(B)

    print()
    print(
        f"Using E=lcm(1,...,{B})"
    )

    print()
    print(
        "p       q       n bits    base    cubic gcd"
    )

    line()

    successes = 0

    for pp, qq in targets:

        nn = pp * qq

        target_hit = False

        for a in range(2, 21):

            if gcd(a,nn) != 1:
                continue

            x = pow(
                a,
                E,
                nn
            )

            g = gcd(
                nn,
                phi3(x)
            )

            if 1 < g < nn:

                successes += 1
                target_hit = True

                print(
                    f"{pp:7d}"
                    f"{qq:8d}"
                    f"{bitlen(nn):10d}"
                    f"{a:8d}"
                    f"{g:18d}"
                )

                break

        if not target_hit:

            print(
                f"{pp:7d}"
                f"{qq:8d}"
                f"{bitlen(nn):10d}"
                f"{'-':>8s}"
                f"{'NONE':>18s}"
            )

    print()
    print(
        f"Successful targets = {successes}/{len(targets)}"
    )


# =============================================================================
# SECTION 11
# NEW QUESTION: CAN WE GENERATE THE ROOT WITH A POWER?
# =============================================================================

def root_generation_equation():

    section(
        "11. ALGEBRAIC POWER-ROOT CONDITION"
    )

    print()
    print(
        "We need:"
    )

    print(
        "    a^(3E) == 1 (mod p)"
    )

    print(
        "but"
    )

    print(
        "    a^E != 1 (mod p)."
    )

    print()
    print(
        "Equivalently, the multiplicative order of a modulo p"
    )

    print(
        "must contain exactly enough 3-adic structure."
    )

    print()
    print(
        "This makes the hidden-factor problem an order-detection problem."
    )

    print()
    print(
        "The experiment therefore records p-1 and q-1 factorization"
    )

    print(
        "for the representative case."
    )

    print()
    print(
        f"p-1 factors = {sp.factorint(P-1)}"
    )

    print(
        f"q-1 factors = {sp.factorint(Q-1)}"
    )


# =============================================================================
# FINAL
# =============================================================================

def final_summary():

    section(
        "12. FINAL RESEARCH DIRECTION"
    )

    print()
    print(
        "The previous polynomial approach is ruled out by:"
    )

    print()
    print(
        "    G(n) == G(0) mod n."
    )

    print()
    print(
        "Therefore the next meaningful class of constructions is:"
    )

    print()
    print(
        "    x = a^E mod n"
    )

    print(
        "    gcd(n, Phi3(x))"
    )

    print()
    print(
        "This may reveal a hidden factor if x has order 3"
    )

    print(
        "modulo one factor but not the other."
    )

    print()
    print(
        "The central experimental question is now:"
    )

    print()
    print(
        "Can cubic-root structure be reached by modular"
    )

    print(
        "exponentiation substantially earlier than ordinary"
    )

    print(
        "factor-order methods such as Pollard p-1?"
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
        "CUBIC ROOT GENERATION BY MODULAR EXPONENTIATION"
    )
    print("=" * 78)

    test_polynomial_obstruction()

    factor_order_data()

    direct_power_search()

    order_control()

    exponent_3_structure()

    smooth_search()

    pollard_p1_comparison()

    group_order_experiment()

    root_of_unity_check()

    multi_target_scaling()

    root_generation_equation()

    final_summary()

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()

