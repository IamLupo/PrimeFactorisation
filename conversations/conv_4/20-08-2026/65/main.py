#!/usr/bin/env python3

import sympy as sp


print("=" * 78)
print("EXPERIMENT 541 START")
print("=" * 78)
print("F2 TRANSLATION DIFFERENCE -> S CHANNEL -> KAPPA QUADRATIC OPERATOR")
print()


# ============================================================================
# SYMBOLS
# ============================================================================

N, S = sp.symbols("N S")
h, k, z = sp.symbols("h k z")
U, V, c = sp.symbols("U V c", nonzero=True)


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
# 1. MODERN KAPPA F2
# ============================================================================

print()
print("[1] F2 DEFINITION")
print("-" * 78)

F2 = 6 * N - S**2 + S

print("  F2 =")
sp.pprint(F2)


# ============================================================================
# 2. DISTINGUISHED TRANSLATION OF N,S
# ============================================================================

print()
print("[2] DISTINGUISHED KAPPA TRANSLATION")
print("-" * 78)

N_h = N + h*S + h**2
S_h = S + 2*h

print("  T_h(N) =")
sp.pprint(N_h)

print()
print("  T_h(S) =")
sp.pprint(S_h)

cert(
    "N translation law",
    N_h - N,
    h*S + h**2,
)

cert(
    "S translation law",
    S_h - S,
    2*h,
)


# ============================================================================
# 3. F2 TRANSLATION
# ============================================================================

print()
print("[3] TRANSLATED F2")
print("-" * 78)

F2_h = sp.expand(
    F2.subs(
        {
            N: N_h,
            S: S_h,
        },
        simultaneous=True,
    )
)

print("  F2(h) =")
sp.pprint(F2_h)

dF2 = sp.factor(F2_h - F2)

print()
print("  Delta_h F2 =")
sp.pprint(dF2)


# ============================================================================
# 4. EXACT S EXTRACTION
# ============================================================================

print()
print("[4] EXACT S EXTRACTION FROM F2 DIFFERENCE")
print("-" * 78)

S_from_F2 = sp.factor(
    dF2 / (2*h) - h - 1
)

print("  S_candidate = Delta_h(F2)/(2h) - h - 1")
print()
print("  candidate =")
sp.pprint(S_from_F2)

cert(
    "S recovered from F2 translation",
    S_from_F2,
    S,
)


# ============================================================================
# 5. S EXTRACTION WITHOUT DIVISION BY h
# ============================================================================

print()
print("[5] DIVISION-FREE S CHANNEL")
print("-" * 78)

division_free = sp.expand(
    dF2 - 2*h*(S + h + 1)
)

cert(
    "Delta_h F2 = 2h(S+h+1)",
    division_free,
    0,
)


# ============================================================================
# 6. FIRST-DIFFERENCE OPERATOR
# ============================================================================

print()
print("[6] FINITE-DIFFERENCE OPERATOR")
print("-" * 78)

D1_F2 = sp.factor(
    dF2
)

D2_F2 = sp.factor(
    (
        F2.subs(
            {
                N: N + 2*h*S + 4*h**2,
                S: S + 4*h,
            },
            simultaneous=True,
        )
        - 2*F2_h
        + F2
    )
)

print("  D(F2) =")
sp.pprint(D1_F2)

print()
print("  D²(F2) =")
sp.pprint(D2_F2)

cert(
    "second difference of F2",
    D2_F2,
    4*h**2,
)


# ============================================================================
# 7. NORMALIZED DIFFERENCE
# ============================================================================

print()
print("[7] NORMALIZED DIFFERENCE CHANNEL")
print("-" * 78)

R_h = sp.factor(
    dF2 / 2 - h*(h + 1)
)

print("  R_h = Delta(F2)/2 - h(h+1)")
sp.pprint(R_h)

cert(
    "normalized difference equals h*S",
    R_h,
    h*S,
)


# ============================================================================
# 8. DIRECT CONSTRUCTION OF THE KAPPA ORBIT
# ============================================================================

print()
print("[8] RECONSTRUCT KAPPA ORBIT FROM F2 DIFFERENCE")
print("-" * 78)

S_extracted = sp.factor(
    R_h / h
)

L_reconstructed = sp.expand(
    N + h*S_extracted + h**2
)

L_target = N + h*S + h**2

print("  S_extracted =")
sp.pprint(S_extracted)

print()
print("  reconstructed L(h) =")
sp.pprint(L_reconstructed)

cert(
    "reconstructed translation quadratic",
    L_reconstructed,
    L_target,
)


# ============================================================================
# 9. SEARCH FOR A UNIVERSAL F2 -> L OPERATOR
# ============================================================================

print()
print("[9] UNIVERSAL OPERATOR ANSATZ")
print("-" * 78)

#
# Search:
#
#   L(h) =
#       a(h)*F2(h)
#     + b(h)*F2(0)
#     + c(h)*N
#     + d(h)*S
#     + e(h)
#
# but we specifically want the operator to recover L without
# explicitly inserting S.
#

a, b, cc, d, e = sp.symbols("a b cc d e")

# Solve for coefficients assumed polynomial/rational in h
# by restricting to the simplest affine-in-h ansatz.

a0, a1 = sp.symbols("a0 a1")
b0, b1 = sp.symbols("b0 b1")
c0, c1 = sp.symbols("c0 c1")
d0, d1 = sp.symbols("d0 d1")
e0, e1, e2 = sp.symbols("e0 e1 e2")

a_h = a0 + a1*h
b_h = b0 + b1*h
c_h = c0 + c1*h
d_h = d0 + d1*h
e_h = e0 + e1*h + e2*h**2

operator_candidate = sp.expand(
    a_h * F2_h
    + b_h * F2
    + c_h * N
    + d_h * S
    + e_h
)

target = N + h*S + h**2

difference = sp.Poly(
    sp.expand(operator_candidate - target),
    N,
    S,
    h,
)

equations = [
    coeff
    for coeff in difference.coeffs()
]

solution = sp.solve(
    equations,
    [
        a0, a1,
        b0, b1,
        c0, c1,
        d0, d1,
        e0, e1, e2,
    ],
    dict=True,
)

print("  affine coefficient solutions =")
sp.pprint(solution)


# ============================================================================
# 10. CHECK WHETHER N-ONLY TERMS ARE ACTUALLY NECESSARY
# ============================================================================

print()
print("[10] F2-ONLY OPERATOR TEST")
print("-" * 78)

a0f, a1f = sp.symbols("a0f a1f")
b0f, b1f = sp.symbols("b0f b1f")
e0f, e1f, e2f = sp.symbols("e0f e1f e2f")

f2_only = sp.expand(
    (a0f + a1f*h) * F2_h
    + (b0f + b1f*h) * F2
    + e0f + e1f*h + e2f*h**2
)

f2_only_poly = sp.Poly(
    sp.expand(f2_only - target),
    N,
    S,
    h,
)

f2_only_equations = f2_only_poly.coeffs()

f2_only_solution = sp.solve(
    f2_only_equations,
    [
        a0f, a1f,
        b0f, b1f,
        e0f, e1f, e2f,
    ],
    dict=True,
)

print("  F2-only affine solutions =")
sp.pprint(f2_only_solution)


# ============================================================================
# 11. HISTORICAL COORDINATE REPRESENTATION
# ============================================================================

print()
print("[11] HISTORICAL REPRESENTATION")
print("-" * 78)

N_hist = (U**2 - V**2) / c**2
S_hist = 2*U / c

F2_hist = sp.factor(
    F2.subs(
        {
            N: N_hist,
            S: S_hist,
        },
        simultaneous=True,
    )
)

print("  F2(U,V) =")
sp.pprint(F2_hist)

print()
print("  F2(U+h*c,V) - F2(U,V) =")

F2_hist_shift = sp.factor(
    F2_hist.subs(U, U + h*c) - F2_hist
)

sp.pprint(F2_hist_shift)


# ============================================================================
# 12. HISTORICAL S EXTRACTION
# ============================================================================

print()
print("[12] HISTORICAL S EXTRACTION")
print("-" * 78)

S_hist_candidate = sp.factor(
    F2_hist_shift / (2*h) - h - 1
)

print("  extracted =")
sp.pprint(S_hist_candidate)

cert(
    "historical S extraction",
    S_hist_candidate,
    S_hist,
)


# ============================================================================
# 13. HISTORICAL N RECONSTRUCTION
# ============================================================================

print()
print("[13] HISTORICAL N RECONSTRUCTION")
print("-" * 78)

N_from_hist = sp.factor(
    F2_hist
    - 2*N_hist
    + S_hist**2
    - S_hist
)

print("  reconstructed N from F2:")
sp.pprint(N_from_hist)

cert(
    "N reconstructed from F2 after adding S²-S",
    N_from_hist,
    4*N_hist,
)


# ============================================================================
# 14. WHY F2 IS NOT ITSELF THE KAPPA QUADRATIC
# ============================================================================

print()
print("[14] F2 VS KAPPA QUADRATIC")
print("-" * 78)

Q_kappa = z**2 - S*z + N

print("  KAPPA quadratic:")
sp.pprint(Q_kappa)

print()
print("  F2 under spectral substitution S=z:")
F2_z = sp.factor(
    F2.subs(S, z)
)

sp.pprint(F2_z)

print()
print("""
  F2 is not the generating quadratic.

  Instead:

      F2 = 6N - S² + S

  is a quadratic functional of the same symmetric coordinate.

  Therefore the operator difference is more informative than F2
  itself.
""")


# ============================================================================
# 15. DISCRIMINANT OF F2 ORBIT
# ============================================================================

print()
print("[15] F2 ORBIT DISCRIMINANT")
print("-" * 78)

F2_poly_h = sp.Poly(
    sp.expand(F2_h),
    h,
)

f2_c2 = sp.factor(
    F2_poly_h.coeff_monomial(h**2)
)

f2_c1 = sp.factor(
    F2_poly_h.coeff_monomial(h)
)

f2_c0 = sp.factor(
    F2_poly_h.coeff_monomial(1)
)

disc_F2_orbit = sp.factor(
    f2_c1**2 - 4*f2_c2*f2_c0
)

print("  h² coefficient =")
sp.pprint(f2_c2)

print()
print("  h coefficient =")
sp.pprint(f2_c1)

print()
print("  constant =")
sp.pprint(f2_c0)

print()
print("  discriminant =")
sp.pprint(disc_F2_orbit)

print()
print("""
  The discriminant of the F2 orbit is NOT expected to equal
  Delta directly.

  The important question is whether it contains a universal
  affine transformation of Delta, because F2 has an extra
  quadratic functional dependence on S.
""")


# ============================================================================
# 16. SOLVE F2 ORBIT DISCRIMINANT RELATION
# ============================================================================

print()
print("[16] DISCRIMINANT RELATION TO DELTA")
print("-" * 78)

Delta = S**2 - 4*N

disc_difference = sp.factor(
    disc_F2_orbit - 16*Delta
)

print("  disc(F2 orbit) - 16*Delta =")
sp.pprint(disc_difference)

# Also test a broader affine relation:
aa, bb = sp.symbols("aa bb")

affine_disc = sp.expand(
    disc_F2_orbit - (aa*Delta + bb)
)

disc_coeffs = sp.Poly(
    affine_disc,
    N,
    S,
).coeffs()

disc_solution = sp.solve(
    disc_coeffs,
    [aa, bb],
    dict=True,
)

print()
print("  affine Delta relation solutions =")
sp.pprint(disc_solution)


# ============================================================================
# 17. OPERATOR COMMUTATOR TEST
# ============================================================================

print()
print("[17] TRANSLATION OPERATOR COMMUTATOR")
print("-" * 78)

def T(expr, amount):
    return sp.expand(
        expr.subs(
            {
                N: N + amount*S + amount**2,
                S: S + 2*amount,
            },
            simultaneous=True,
        )
    )


commutator = sp.factor(
    T(T(F2, h), k) - T(T(F2, k), h)
)

print("  T_h T_k(F2) - T_k T_h(F2) =")
sp.pprint(commutator)

cert(
    "translation operators commute on F2",
    commutator,
    0,
)


# ============================================================================
# 18. GROUP LAW
# ============================================================================

print()
print("[18] TRANSLATION GROUP LAW")
print("-" * 78)

group_difference = sp.factor(
    T(F2, h).subs(
        {
            N: N + k*S + k**2,
            S: S + 2*k,
        },
        simultaneous=True,
    )
    - T(F2, h + k)
)

print("  T_k(T_h(F2)) - T_(h+k)(F2) =")
sp.pprint(group_difference)

cert(
    "translation group law",
    group_difference,
    0,
)


# ============================================================================
# 19. FINAL STRUCTURAL CERTIFICATE
# ============================================================================

print()
print("[19] STRUCTURAL CERTIFICATE")
print("-" * 78)

print("""
The exact F2 translation law is:

    F2 = 6N - S² + S

    T_h(N) = N + hS + h²
    T_h(S) = S + 2h

therefore:

    T_h(F2)-F2
        = 2h(S+h+1).

Hence:

    S
      = [T_h(F2)-F2]/(2h) - h - 1.

This is a direct operator extraction of the symmetric channel.

The important distinction is:

    F2 itself
        does not equal
    N + hS + h².

Instead, the FIRST DIFFERENCE of F2 contains S.

Therefore the upstream bridge has potentially changed form:

    homogeneous object F2
          |
          v
    translation difference
          |
          v
          S
          |
          +---- N
          |
          v
    Q(z)=z²-Sz+N.

The decisive unresolved question is no longer whether F2 has
the KAPPA quadratic orbit.

It does not.

The question is whether the ORIGINAL homogeneous-layer construction
provides an intrinsic operator T for which this first difference
is available without already knowing S.

If yes, F2 itself may be the upstream S-channel carrier.
""")


# ============================================================================
# 20. FINAL AUDIT
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 541 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
The next experiment should NOT introduce another synthetic
quadratic form.

Instead, take the actual homogeneous-layer index/operator and
test whether it realizes the transformation:

    (N,S)
        ->
    (N+hS+h², S+2h).

Then ask whether the already-existing F2 layer transforms by:

    T_h(F2)-F2 = 2h(S+h+1).

If the operator exists intrinsically upstream, extract:

    S = Delta(F2)/(2h)-h-1

and immediately construct:

    Q(z)=z²-Sz+N.

The important target is therefore the OPERATOR that produces the
F2 layer shift, not another expression for F2 itself.

Do not use factor pairs.
Do not use continued fractions.
Do not fit numerical data.
Do not insert p or q.
""")
