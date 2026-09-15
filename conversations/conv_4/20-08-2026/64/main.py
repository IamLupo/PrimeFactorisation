#!/usr/bin/env python3

import sympy as sp


print("=" * 78)
print("EXPERIMENT 540 START")
print("=" * 78)
print("KAPPA HOMOGENEOUS LAYER -> HISTORICAL CONIC TRANSLATION CHANNEL")
print()


# ============================================================================
# SYMBOLS
# ============================================================================

S, N = sp.symbols("S N")
h = sp.symbols("h")
U, V, c = sp.symbols("U V c", nonzero=True)

x, y = sp.symbols("x y")


# ============================================================================
# HELPERS
# ============================================================================

failures = 0


def simp(expr):
    return sp.factor(sp.expand(sp.simplify(expr)))


def cert(label, expr, expected=0):
    global failures

    diff = simp(expr - expected)

    ok = diff == 0

    print(f"  {label}")
    print(f"    difference = {diff}")
    print(f"    PASS = {ok}")

    if not ok:
        failures += 1

    return ok


# ============================================================================
# 1. MODERN KAPPA LAYERS
# ============================================================================

print("[1] MODERN KAPPA HOMOGENEOUS LAYERS")
print("-" * 78)

# Established KAPPA quantities.
#
# F2 = 6N - S^2 + S
# F3 = (S+1) F2
#
# The experiment treats these as algebraic layer objects.

F2 = 6 * N - S**2 + S
F3 = (S + 1) * F2

print("  F2 =")
sp.pprint(F2)

print()
print("  F3 =")
sp.pprint(F3)

cert(
    "F3=(S+1)F2",
    F3,
    (S + 1) * F2,
)


# ============================================================================
# 2. SUBSTITUTE HISTORICAL COORDINATES
# ============================================================================

print()
print("[2] HISTORICAL COORDINATE SUBSTITUTION")
print("-" * 78)

# Historical normalized coordinates:
#
#     S = 2U/c
#     N = (U^2 - V^2)/c^2

S_hist = 2 * U / c
N_hist = (U**2 - V**2) / c**2

F2_hist = sp.factor(
    F2.subs(
        {
            S: S_hist,
            N: N_hist,
        },
        simultaneous=True,
    )
)

F3_hist = sp.factor(
    F3.subs(
        {
            S: S_hist,
            N: N_hist,
        },
        simultaneous=True,
    )
)

print("  F2(U,V) =")
sp.pprint(F2_hist)

print()
print("  F3(U,V) =")
sp.pprint(F3_hist)


# ============================================================================
# 3. U-TRANSLATION OPERATOR
# ============================================================================

print()
print("[3] DISTINGUISHED U-TRANSLATION")
print("-" * 78)

def T_U(expr, amount=1):
    return sp.expand(
        expr.subs(U, U + amount * c)
    )


F2_shift = sp.factor(T_U(F2_hist, h))
F3_shift = sp.factor(T_U(F3_hist, h))

print("  T_h(U)=U+h*c")
print()
print("  T_h(F2) =")
sp.pprint(F2_shift)

print()
print("  T_h(F3) =")
sp.pprint(F3_shift)


# ============================================================================
# 4. FIRST TRANSLATION RESPONSE
# ============================================================================

print()
print("[4] FIRST TRANSLATION RESPONSE")
print("-" * 78)

dF2 = sp.factor(
    T_U(F2_hist, h) - F2_hist
)

dF3 = sp.factor(
    T_U(F3_hist, h) - F3_hist
)

print("  Delta_h F2 =")
sp.pprint(dF2)

print()
print("  Delta_h F3 =")
sp.pprint(dF3)


# Expected first response from S -> S+2h:
#
#     F2(S+2h,N+hS+h^2) - F2(S,N)
#
# but this must be derived symbolically.

S_shift = S + 2 * h
N_shift = N + h * S + h**2

F2_SN_shift = sp.expand(
    F2.subs(
        {
            S: S_shift,
            N: N_shift,
        },
        simultaneous=True,
    )
)

expected_dF2_SN = sp.factor(
    F2_SN_shift - F2
)

print()
print("  F2 response from N,S translation =")
sp.pprint(expected_dF2_SN)

cert(
    "historical and N,S translation responses agree",
    dF2,
    expected_dF2_SN.subs(
        {
            S: S_hist,
            N: N_hist,
        },
        simultaneous=True,
    ),
)


# ============================================================================
# 5. SECOND DIFFERENCE OF F2
# ============================================================================

print()
print("[5] SECOND DIFFERENCE OF F2")
print("-" * 78)

F2_0 = F2_hist
F2_1 = T_U(F2_hist, 1)
F2_2 = T_U(F2_hist, 2)

second_F2 = sp.factor(
    F2_2 - 2 * F2_1 + F2_0
)

print("  Delta^2 F2 =")
sp.pprint(second_F2)


# ============================================================================
# 6. COMPARE WITH HISTORICAL v TRANSLATION
# ============================================================================

print()
print("[6] HISTORICAL v-KERNEL TRANSLATION")
print("-" * 78)

v = y**2 - x**2 + 3 * x - 2

# Historical coordinates:
#
#     U = 2y
#     V = 2x - 3

v_hist = sp.expand(
    v.subs(
        {
            x: (V + 3) / 2,
            y: U / 2,
        },
        simultaneous=True,
    )
)

print("  v(U,V) =")
sp.pprint(v_hist)

cert(
    "4v-1 = U^2-V^2",
    4 * v_hist - 1,
    U**2 - V**2,
)


# ============================================================================
# 7. TRANSLATED HISTORICAL v
# ============================================================================

print()
print("[7] TRANSLATED HISTORICAL v")
print("-" * 78)

v_shift = sp.factor(
    v_hist.subs(U, U + h * c) - v_hist
)

print("  T_h(v)-v =")
sp.pprint(v_shift)

expected_v_shift = h * c * U / 2 + h**2 * c**2 / 4

cert(
    "historical v translation",
    v_shift,
    expected_v_shift,
)


# ============================================================================
# 8. F2 IN TERMS OF v
# ============================================================================

print()
print("[8] F2 -> HISTORICAL CONIC KERNEL")
print("-" * 78)

K = U**2 - V**2

F2_K = sp.factor(
    F2_hist
)

print("  K = U^2-V^2")
print()
print("  F2(K,U) =")
sp.pprint(F2_K)


# Test whether F2 can be decomposed into powers of K and U.
#
# This is important because K supplies N while U supplies S.

print()
print("  F2 expanded in U,V:")
sp.pprint(sp.expand(F2_K))


# ============================================================================
# 9. TRANSLATION RESPONSE DECOMPOSITION
# ============================================================================

print()
print("[9] TRANSLATION RESPONSE DECOMPOSITION")
print("-" * 78)

response = sp.expand(
    dF2
)

response_poly = sp.Poly(response, h)

resp_h2 = sp.factor(
    response_poly.coeff_monomial(h**2)
)

resp_h1 = sp.factor(
    response_poly.coeff_monomial(h)
)

print("  coefficient of h =")
sp.pprint(resp_h1)

print()
print("  coefficient of h^2 =")
sp.pprint(resp_h2)


# ============================================================================
# 10. S-CHANNEL TEST
# ============================================================================

print()
print("[10] PURE S-CHANNEL TEST")
print("-" * 78)

# The historical translation has:
#
#     S = 2U/c.
#
# We ask whether the h-linear response contains an affine
# multiple of S and whether the remaining part depends only
# on N/invariant quantities.

linear_response = resp_h1

print("  linear translation channel =")
sp.pprint(linear_response)

S_channel = sp.factor(
    linear_response.subs(U, S * c / 2)
)

print()
print("  after U=(S*c)/2:")
sp.pprint(S_channel)


# ============================================================================
# 11. FACTOR OUT THE HISTORICAL MOVING COORDINATE
# ============================================================================

print()
print("[11] MOVING / INVARIANT SEPARATION")
print("-" * 78)

# Write the linear response as:
#
#     A(V,c)*U + B(V,c)
#
# and inspect the U coefficient.

linear_U_coeff = sp.factor(
    sp.diff(linear_response, U)
)

linear_residual = sp.factor(
    linear_response - linear_U_coeff * U
)

print("  U coefficient =")
sp.pprint(linear_U_coeff)

print()
print("  U-independent residual =")
sp.pprint(linear_residual)


# ============================================================================
# 12. DISCRIMINANT RECONSTRUCTION FROM THE F2 ORBIT
# ============================================================================

print()
print("[12] F2 TRANSLATION ORBIT AS A POLYNOMIAL IN h")
print("-" * 78)

F2_orbit_poly = sp.Poly(
    sp.expand(F2_shift),
    h,
)

orbit_a2 = sp.factor(
    F2_orbit_poly.coeff_monomial(h**2)
)

orbit_a1 = sp.factor(
    F2_orbit_poly.coeff_monomial(h)
)

orbit_a0 = sp.factor(
    F2_orbit_poly.coeff_monomial(1)
)

print("  coefficient h^2 =")
sp.pprint(orbit_a2)

print()
print("  coefficient h =")
sp.pprint(orbit_a1)

print()
print("  constant =")
sp.pprint(orbit_a0)


# ============================================================================
# 13. ORBIT DISCRIMINANT
# ============================================================================

print()
print("[13] F2-ORBIT DISCRIMINANT")
print("-" * 78)

if orbit_a2 != 0:
    orbit_disc = sp.factor(
        orbit_a1**2 - 4 * orbit_a2 * orbit_a0
    )
else:
    orbit_disc = sp.Integer(0)

print("  discriminant in h =")
sp.pprint(orbit_disc)


# Compare with the historical gap channel.
target_gap = sp.factor(
    4 * V**2 / c**2
)

print()
print("  target historical gap =")
sp.pprint(target_gap)


# ============================================================================
# 14. NORMALIZED F2 ORBIT
# ============================================================================

print()
print("[14] NORMALIZED ORBIT RECONSTRUCTION")
print("-" * 78)

# F2 itself is not the KAPPA quadratic.
#
# Search for an affine normalization:
#
#     G_h = a*F2_h + b
#
# whose h^2 coefficient is 1.
#
# Then inspect its h coefficient and discriminant.

scale = sp.symbols("scale", nonzero=True)

G = sp.expand(
    scale * F2_shift
)

G_poly = sp.Poly(G, h)

G_h2 = sp.factor(
    G_poly.coeff_monomial(h**2)
)

G_h1 = sp.factor(
    G_poly.coeff_monomial(h)
)

print("  scaled h^2 coefficient =")
sp.pprint(G_h2)

print()
print("  scaled h coefficient =")
sp.pprint(G_h1)


# Determine scale needed to normalize h^2 to one.
scale_solution = sp.solve(
    sp.Eq(G_h2, 1),
    scale,
)

print()
print("  scale values giving monic orbit =")
sp.pprint(scale_solution)


# ============================================================================
# 15. DIRECT COMPARISON WITH KAPPA ORBIT
# ============================================================================

print()
print("[15] DIRECT COMPARISON WITH KAPPA ORBIT")
print("-" * 78)

L_kappa = (
    h**2
    + S_hist * h
    + N_hist
)

print("  KAPPA orbit L(h) =")
sp.pprint(L_kappa)

# Try the simplest canonical KAPPA expression directly.
Kappa_from_hist = sp.expand(L_kappa)

# Verify its translation form.
cert(
    "historical KAPPA orbit",
    Kappa_from_hist,
    (U + h*c)**2 / c**2 - V**2 / c**2,
)


# ============================================================================
# 16. F2 DIFFERENCE FROM THE KAPPA ORBIT
# ============================================================================

print()
print("[16] DIFFERENCE BETWEEN F2 AND TRANSLATION KERNEL")
print("-" * 78)

difference_F2_L = sp.factor(
    F2_hist - L_kappa
)

print("  F2 - L(h=0) =")
sp.pprint(difference_F2_L)

print()
print("""
  This is the important structural test:

      Is F2 merely a function of the translation quadratic L?

  If it is affine in L, then the homogeneous KAPPA layer is
  directly carrying the historical translation orbit.

  If not, F2 belongs to a higher layer and the useful bridge
  may require a specific combination of layers.
""")


# ============================================================================
# 17. SEARCH LOW-DEGREE COMBINATIONS OF F2, F3 AND L
# ============================================================================

print()
print("[17] LOW-DEGREE OPERATOR COMBINATIONS")
print("-" * 78)

# Search simple expressions:
#
#     F2 + a*N + b*S + c
#     F3 + a*F2 + b*N + c*S + d
#
# for coefficients that become especially simple under
# historical translation.

a, b, cc, d = sp.symbols("a b cc d")

candidate1 = sp.expand(
    F2 + a*N + b*S + cc
)

candidate1_hist = sp.expand(
    candidate1.subs(
        {
            S: S_hist,
            N: N_hist,
        },
        simultaneous=True,
    )
)

print("  candidate family:")
sp.pprint(candidate1_hist)

candidate2 = sp.expand(
    F3 + a*F2 + b*N + cc*S + d
)

candidate2_hist = sp.expand(
    candidate2.subs(
        {
            S: S_hist,
            N: N_hist,
        },
        simultaneous=True,
    )
)

print()
print("  second candidate family:")
sp.pprint(candidate2_hist)


# ============================================================================
# 18. TRANSLATION OF N AND S
# ============================================================================

print()
print("[18] EXACT N,S TRANSLATION")
print("-" * 78)

N_h = sp.expand(
    N_hist.subs(U, U + h*c)
)

S_h = sp.expand(
    S_hist.subs(U, U + h*c)
)

print("  T_h(N) =")
sp.pprint(N_h)

print()
print("  T_h(S) =")
sp.pprint(S_h)

cert(
    "N translation",
    N_h - N_hist,
    h*S_hist + h**2,
)

cert(
    "S translation",
    S_h - S_hist,
    2*h,
)


# ============================================================================
# 19. TRANSLATION OPERATOR ON F2
# ============================================================================

print()
print("[19] OPERATOR ACTION IN (N,S)")
print("-" * 78)

F2_transformed_NS = sp.expand(
    F2.subs(
        {
            N: N + h*S + h**2,
            S: S + 2*h,
        },
        simultaneous=True,
    )
)

print("  T_h(F2) =")
sp.pprint(F2_transformed_NS)

cert(
    "F2 translation from N,S laws",
    F2_transformed_NS,
    F2_SN_shift,
)


# ============================================================================
# 20. DIFFERENCE OPERATOR ON HOMOGENEOUS LAYER
# ============================================================================

print()
print("[20] HOMOGENEOUS DIFFERENCE OPERATOR")
print("-" * 78)

D_F2 = sp.factor(
    F2_transformed_NS - F2
)

D2_F2 = sp.factor(
    F2_transformed_NS.subs(
        {
            N: N + h*S + h**2,
            S: S + 2*h,
        },
        simultaneous=True,
    )
    - 2*F2_transformed_NS
    + F2
)

print("  D(F2) =")
sp.pprint(D_F2)

print()
print("  D²(F2) =")
sp.pprint(D2_F2)


# ============================================================================
# 21. FINAL STRUCTURAL INTERPRETATION
# ============================================================================

print()
print("[21] STRUCTURAL INTERPRETATION")
print("-" * 78)

print("""
The previous experiments established that an arbitrary quadratic
form under translation automatically has:

    quadratic orbit
    constant second difference
    linear first difference.

Experiment 540 asks a sharper question about the REAL KAPPA layer:

    F2 = 6N - S² + S.

Under the distinguished translation:

    S -> S + 2h

    N -> N + hS + h².

Therefore the relevant question is whether F2 itself, or a
low-complexity combination of the existing homogeneous layers,
contains the same translation quadratic

    L(h) = N + hS + h².

The strongest possible outcome would be an exact operator identity

    L(h) = A(F2,T(F2),...) + invariant

with the invariant part depending only on the transverse channel.

The experiment therefore does NOT search for another factorization
identity.

It tests whether the already-existing KAPPA homogeneous layer
responds to the distinguished historical translation operator in
the same algebraic way as the historical conic.
""")


# ============================================================================
# 22. FINAL AUDIT
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 540 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
If F2 or a low-degree combination of the existing homogeneous
layers produces the translation orbit, the next step is to solve
the resulting operator equation exactly and identify the invariant
remainder.

If it does NOT, that is also informative:

    the historical/KAPPA quadratic orbit is not contained in F2
    alone and must arise from a different homogeneous layer or
    from a specific layer-combination operator.

The next productive search should therefore remain confined to
the already-established KAPPA homogeneous layers and their exact
translation operator.

Do not return to continued fractions.
Do not enumerate factor pairs.
Do not fit numerical data.
Do not insert p or q.
""")
