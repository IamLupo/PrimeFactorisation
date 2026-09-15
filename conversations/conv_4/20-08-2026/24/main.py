#!/usr/bin/env python3

import math
import sympy as sp

print("EXPERIMENT 501 START")
print("=" * 78)
print("SHIFTED-PRODUCT GCD / CONGRUENCE FACTOR-LEAK SEARCH")
print("=" * 78)
print()

# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")
x, a, b, z = sp.symbols("x a b z")


def fact(expr):
    return sp.factor(sp.expand(expr))


def zero(expr):
    return sp.expand(sp.cancel(expr)) == 0


# ============================================================================
# CORE DEFINITIONS
# ============================================================================

M_pq = (p + x) * (q + x)
M_NS = N + S*x + x**2

Delta = S**2 - 4*N


# ============================================================================
# 1. SHIFTED PRODUCT IDENTITY
# ============================================================================

print("[1] SHIFTED-PRODUCT IDENTITY")
print("-" * 78)

identity = sp.expand(
    M_pq.subs(
        {
            p + q: S,
            p*q: N,
        }
    ) - M_NS
)

# Direct substitution p+q is not reliable inside arbitrary SymPy expressions,
# so verify using the explicit difference first.
direct_difference = sp.expand(
    (p + x)*(q + x) - (p*q + (p+q)*x + x**2)
)

print("  (p+x)(q+x) - pq - (p+q)x - x^2 =")
print("   ", fact(direct_difference))
print("  PASS =", zero(direct_difference))
print()


# ============================================================================
# 2. SHIFTED PRODUCT MODULO N
# ============================================================================

print("[2] SHIFTED PRODUCT MODULO N")
print("-" * 78)

mod_N_difference = sp.expand(
    (p + x)*(q + x) - (p*q + x*(p+q) + x**2)
)

print("  M_x = N + S*x + x^2")
print()
print("  M_x - x^2 = N + S*x")
print()
print("  Therefore modulo N:")
print("    M_x - x^2 == S*x (mod N)")
print()

print(
    "  Symbolic certificate:",
    fact(mod_N_difference),
)
print(
    "  PASS =",
    zero(mod_N_difference),
)
print()


# ============================================================================
# 3. GCD CHANNELS
# ============================================================================

print("[3] SYMBOLIC GCD CHANNELS")
print("-" * 78)

print("For integer x:")
print()
print("  gcd(M_x - x^2, N) = gcd(S*x, N)")
print("  gcd(M_x - x^2, N) contains a factor r of N")
print("  only when r | S*x.")
print()

print("For gcd(x,N)=1 this reduces to:")
print()
print("  gcd(M_x-x^2,N) = gcd(S,N)")
print()

print("Since S=p+q:")
print()
print("  gcd(S,N) = gcd(p+q,pq)")
print()

print(
    "For distinct primes p,q this is normally 1, "
    "so this first channel is expected to be weak."
)
print()


# ============================================================================
# 4. SPECIAL +1 / -1 CHANNELS
# ============================================================================

print("[4] SPECIAL SHIFTS +1 AND -1")
print("-" * 78)

M1 = (p + 1)*(q + 1)
Mm1 = (p - 1)*(q - 1)

M1_NS = N + S + 1
Mm1_NS = N - S + 1

print("  M(+1) =", M1_NS)
print("  M(-1) =", Mm1_NS)
print()

special_product = sp.expand(
    M1_NS * Mm1_NS
)

special_square = sp.expand(
    (N + 1)**2 - special_product
)

special_delta = sp.expand(
    special_square - 4*N
)

print("  (N+1)^2 - M(+1)M(-1) =", fact(special_square))
print("  target S^2 =", fact(S**2))
print(
    "  S^2 difference =",
    fact(special_square - S**2),
)
print(
    "  S^2 PASS =",
    zero(special_square - S**2),
)
print()

print(
    "  Delta = (N+1)^2 - M(+1)M(-1) - 4N"
)
print(
    "  recovered Delta =",
    fact(special_delta),
)
print(
    "  target Delta =",
    fact(Delta),
)
print(
    "  Delta difference =",
    fact(special_delta - Delta),
)
print(
    "  Delta PASS =",
    zero(special_delta - Delta),
)
print()


# ============================================================================
# 5. GCD SEARCH OVER NATURAL SHIFTS
# ============================================================================

print("[5] NUMERICAL GCD SEARCH OVER NATURAL SHIFTS")
print("-" * 78)

instances = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),
]

shifts = [
    -20, -16, -12, -10, -8, -6, -5, -4, -3, -2, -1,
    1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20
]

interesting_hits = []

for pp, qq in instances:
    nn = pp * qq
    ss = pp + qq

    print(f"  ({pp},{qq}) N={nn}")

    local_hits = []

    for xx in shifts:
        mx = (pp + xx) * (qq + xx)

        # Primary congruence channel
        g1 = math.gcd(
            mx - xx**2,
            nn,
        )

        # Variants
        g2 = math.gcd(
            mx - xx**2 - xx,
            nn,
        )

        g3 = math.gcd(
            mx - xx**2 + xx,
            nn,
        )

        # Remove trivial gcd=1 and gcd=N cases.
        for label, g in [
            ("M-x^2", g1),
            ("M-x^2-x", g2),
            ("M-x^2+x", g3),
        ]:
            if 1 < g < nn:
                local_hits.append(
                    (xx, label, g)
                )

    if local_hits:
        for hit in local_hits:
            print(
                f"    HIT x={hit[0]:>3} "
                f"{hit[1]:>12} gcd={hit[2]}"
            )
            interesting_hits.append(
                (pp, qq, hit)
            )
    else:
        print("    no nontrivial gcd hits")

print()


# ============================================================================
# 6. CHECK THE PRIMARY CHANNEL ANALYTICALLY
# ============================================================================

print("[6] PRIMARY GCD CHANNEL ANALYTIC CHECK")
print("-" * 78)

print("For any integer x:")
print()
print("  M_x - x^2 = N + S*x")
print()
print("Hence:")
print()
print("  gcd(M_x-x^2,N) = gcd(S*x,N)")
print()

print("For gcd(x,N)=1:")
print()
print("  gcd(M_x-x^2,N) = gcd(S,N)")
print()

print(
    "This means the natural shifted-product congruence does not "
    "normally expose p or q directly."
)
print()


# ============================================================================
# 7. SEARCH FOR SHIFTS THAT FORCE A FACTOR
# ============================================================================

print("[7] SHIFT CONDITIONS FOR FACTOR EXPOSURE")
print("-" * 78)

print("Modulo p:")
print()
print("  M_x = (p+x)(q+x)")
print("      == x(q+x) (mod p)")
print()
print("Therefore p divides M_x exactly when:")
print()
print("  x == 0 (mod p)")
print("  or")
print("  x == -q (mod p)")
print()

print("Modulo q:")
print()
print("  M_x == x(p+x) (mod q)")
print()
print("Therefore q divides M_x exactly when:")
print()
print("  x == 0 (mod q)")
print("  or")
print("  x == -p (mod q)")
print()

print(
    "The nontrivial conditions require knowledge related to the "
    "unknown complementary factor."
)
print()


# ============================================================================
# 8. TWO-SHIFT GCD CHANNEL
# ============================================================================

print("[8] TWO-SHIFT DIFFERENCE CHANNEL")
print("-" * 78)

M_a = (p + a)*(q + a)
M_b = (p + b)*(q + b)

shift_difference = sp.expand(
    M_a - M_b
)

expected_difference = sp.expand(
    (a - b)*(p + q + a + b)
)

print(
    "  M_a - M_b =",
    fact(shift_difference),
)
print(
    "  expected   =",
    fact(expected_difference),
)
print(
    "  difference =",
    fact(shift_difference - expected_difference),
)
print(
    "  PASS =",
    zero(shift_difference - expected_difference),
)
print()


# ============================================================================
# 9. TWO-SHIFT MOD N CHANNEL
# ============================================================================

print("[9] TWO-SHIFT MOD-N CHANNEL")
print("-" * 78)

print("Modulo N:")
print()
print("  M_a == a(S+a)")
print("  M_b == b(S+b)")
print()

print("Subtracting:")
print()
print("  M_a-M_b")
print("      == (a-b)S + (a^2-b^2)")
print()
print("      == (a-b)(S+a+b)")
print()

print(
    "Thus, when gcd(a-b,N)=1, the difference channel "
    "contains S modulo N."
)
print()


# ============================================================================
# 10. NUMERICAL TWO-SHIFT GCD SEARCH
# ============================================================================

print("[10] TWO-SHIFT GCD SEARCH")
print("-" * 78)

pair_shifts = [
    (-3, -2),
    (-2, -1),
    (-1, 1),
    (1, 2),
    (2, 3),
    (1, 3),
    (-3, 3),
]

pair_hits = []

for pp, qq in instances:
    nn = pp * qq

    print(f"  ({pp},{qq})")

    local = []

    for aa, bb in pair_shifts:
        ma = (pp + aa)*(qq + aa)
        mb = (pp + bb)*(qq + bb)

        d = ma - mb

        candidates = {
            "diff": d,
            "diff-corrected":
                d - (aa*aa - bb*bb),
            "diff-div": d // (aa - bb)
            if (aa - bb) != 0 else d,
        }

        for label, value in candidates.items():
            g = math.gcd(value, nn)

            if 1 < g < nn:
                local.append(
                    (aa, bb, label, g)
                )

    if local:
        for hit in local:
            print(
                f"    HIT ({hit[0]},{hit[1]}) "
                f"{hit[2]} gcd={hit[3]}"
            )
            pair_hits.append(
                (pp, qq, hit)
            )
    else:
        print("    no nontrivial pair-gcd hits")

print()


# ============================================================================
# 11. DISCRIMINANT AS ROOT-GAP
# ============================================================================

print("[11] ROOT-GAP / DISCRIMINANT CHANNEL")
print("-" * 78)

quadratic = z**2 - S*z + N

disc_quadratic = sp.discriminant(
    quadratic,
    z,
)

print("  quadratic =", fact(quadratic))
print("  discriminant =", fact(disc_quadratic))
print("  Delta =", fact(Delta))
print(
    "  difference =",
    fact(disc_quadratic - Delta),
)
print(
    "  PASS =",
    zero(disc_quadratic - Delta),
)
print()

print(
    "  roots = (S +/- sqrt(Delta))/2"
)
print(
    "  factor gap = sqrt(Delta) = |p-q|"
)
print()


# ============================================================================
# 12. NUMERICAL DISCRIMINANT AUDIT
# ============================================================================

print("[12] NUMERICAL DISCRIMINANT AUDIT")
print("-" * 78)

disc_failures = 0

for pp, qq in instances:
    nn = pp * qq
    ss = pp + qq

    d = ss**2 - 4*nn
    expected = (pp - qq)**2

    passed = d == expected

    if not passed:
        disc_failures += 1

    print(
        f"  ({pp},{qq}) "
        f"Delta={d} "
        f"expected={(pp-qq)**2} "
        f"PASS={passed}"
    )

print()
print(
    "  DISCRIMINANT FAILURES =",
    disc_failures,
)
print()


# ============================================================================
# 13. INFORMATION-MODEL BOUNDARY
# ============================================================================

print("[13] INFORMATION-MODEL BOUNDARY")
print("-" * 78)

print("The exact downstream maps are:")
print()
print("  M_x")
print("    -> S = (M_x-N-x^2)/x")
print("    -> z^2-Sz+N")
print("    -> p,q")
print()

print("With opposite shifts:")
print()
print("  M(+x), M(-x)")
print("    -> S^2")
print("    -> Delta=S^2-4N")
print("    -> S")
print("    -> p,q")
print()

print("The unresolved question remains:")
print()
print("  N-only construction")
print("      -> M_x / Delta / gap-square")
print("      -> S")
print("      -> p,q")
print()


# ============================================================================
# 14. RESULTS
# ============================================================================

print("[14] EXPERIMENT STATUS")
print("-" * 78)

identity_pass = zero(direct_difference)

s2_special_pass = zero(
    special_square - S**2
)

delta_special_pass = zero(
    special_delta - Delta
)

primary_disc_pass = zero(
    disc_quadratic - Delta
)

print(
    "  shifted-product identity:",
    identity_pass,
)

print(
    "  (N+1)^2-M(+1)M(-1)=S^2:",
    s2_special_pass,
)

print(
    "  corrected Delta identity:",
    delta_special_pass,
)

print(
    "  quadratic discriminant:",
    primary_disc_pass,
)

print(
    "  numerical discriminant audit:",
    disc_failures == 0,
)

print(
    "  nontrivial single-shift gcd hits:",
    len(interesting_hits),
)

print(
    "  nontrivial two-shift gcd hits:",
    len(pair_hits),
)

print()

print("MAIN RESULT")
print("-" * 78)

print(
    "The shifted-product algebra now gives two distinct channels:"
)
print()
print("  1. Symmetric channel:")
print("       M(+x)+M(-x) = 2(N+x^2)")
print()
print("  2. Antisymmetric channel:")
print("       M(+x)-M(-x) = 2Sx")
print()

print(
    "The opposite-shift product gives S^2:"
)
print()
print(
    "  S^2 = ((N+x^2)^2-M(+x)M(-x))/x^2"
)
print()

print(
    "and therefore the actual gap-square is:"
)
print()
print(
    "  Delta=(p-q)^2"
)
print(
    "       = ((N+x^2)^2-M(+x)M(-x))/x^2 - 4N"
)
print()

print(
    "For x=1:"
)
print(
    "  S^2 = (N+1)^2-M(+1)M(-1)"
)
print(
    "  Delta = (N+1)^2-M(+1)M(-1)-4N"
)
print()

print(
    "The GCD search tests whether this shifted-product information "
    "leaks a factor directly through modular structure."
)
print()

print("NEXT TARGET")
print("-" * 78)

print(
    "If the GCD channels produce no nontrivial factors, "
    "the next upstream experiment should stop using M_x entirely "
    "and inspect the existing N-only homogeneous-layer quantities "
    "for an object whose square or resultant has the signature"
)
print()
print("  S^2 - 4N")
print()
print(
    "rather than trying to reconstruct S first."
)

print()
print("=" * 78)
print("EXPERIMENT 501 FINISHED")
print("=" * 78)
