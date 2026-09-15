#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 564
# ==============================================================================
# F2 AS A QUADRATIC IN THE 2020 ROOT COORDINATES
#
# Goal:
#
#   The 2020 construction already gives
#
#       P = A/r
#       Q = B/r
#
#   and therefore
#
#       N     = P Q
#       S     = P + Q
#       Delta = (P-Q)^2.
#
#   The known homogeneous layer is
#
#       F2 = 4 P Q - P^2 - Q^2 + P + Q.
#
#   Instead of trying to recover S from N alone or inventing F4/F5/F6,
#   treat F2 itself as a quadratic observable in the ROOT coordinates.
#
#   The experiment asks:
#
#       1. What quadratic equation does F2 impose on P?
#       2. What is its exact discriminant?
#       3. Is the discriminant expressible through Q/N/S/Delta?
#       4. Does solving the F2 equation preserve the 2020 root structure?
#       5. Can the same construction be written directly in x,y?
#
#   No numerical fitting.
#   No factor enumeration.
#   No F4/F5/F6 assumptions.
# ==============================================================================


print("=" * 78)
print("EXPERIMENT 564 START")
print("=" * 78)
print("F2 AS A QUADRATIC IN THE 2020 ROOT COORDINATES")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

A, B, r = sp.symbols("A B r", nonzero=True)
P, Q = sp.symbols("P Q")
x, y, z = sp.symbols("x y z")

failures = 0


# ==============================================================================
# HELPERS
# ==============================================================================

def clean(expr):
    return sp.factor(sp.cancel(sp.together(expr)))


def numerator(expr):
    return sp.factor(sp.together(expr).as_numer_denom()[0])


def check(name, lhs, rhs=0):
    global failures

    diff = clean(lhs - rhs)
    num = numerator(diff)

    passed = sp.expand(num) == 0

    print(f"  {name}")
    print("    numerator difference =")
    sp.pprint(num)
    print(f"    PASS = {passed}")
    print()

    if not passed:
        failures += 1

    return passed


def show(name, expr):
    print(f"  {name} =")
    sp.pprint(clean(expr))
    print()


# ==============================================================================
# [1] 2020 ROOT COORDINATES
# ==============================================================================

print("[1] 2020 ROOT COORDINATES")
print("-" * 78)

Pxy = clean(A / r)
Qxy = clean(B / r)

N = clean(Pxy * Qxy)
S = clean(Pxy + Qxy)
Delta = clean((Pxy - Qxy)**2)

show("P=A/r", Pxy)
show("Q=B/r", Qxy)
show("N=PQ", N)
show("S=P+Q", S)
show("Delta=(P-Q)^2", Delta)

check("N=PQ", N, Pxy * Qxy)
check("S=P+Q", S, Pxy + Qxy)
check("Delta=(P-Q)^2", Delta, (Pxy - Qxy)**2)

print()


# ==============================================================================
# [2] KNOWN F2 LAYER
# ==============================================================================

print("[2] KNOWN F2 LAYER")
print("-" * 78)

F2_PQ = clean(
    4*P*Q - P**2 - Q**2 + P + Q
)

F2_NSDelta = clean(
    2*N + S - Delta
)

F2_root = clean(
    F2_PQ.subs({
        P: Pxy,
        Q: Qxy,
    })
)

show("F2(P,Q)", F2_PQ)
show("F2=2N+S-Delta", F2_NSDelta)
show("F2(A,B,r)", F2_root)

check(
    "F2(P,Q)=2N+S-Delta",
    F2_PQ.subs({P: Pxy, Q: Qxy}),
    F2_NSDelta,
)

print()


# ==============================================================================
# [3] F2 IN SUM/GAP COORDINATES
# ==============================================================================

print("[3] F2 IN SUM/GAP COORDINATES")
print("-" * 78)

sum_var, gap_var = sp.symbols("sum_var gap_var")

F2_sum_gap = clean(
    sp.expand(
        F2_PQ.subs({
            P: (sum_var + gap_var)/2,
            Q: (sum_var - gap_var)/2,
        })
    )
)

expected_sum_gap = clean(
    (sum_var**2 + 2*sum_var - 3*gap_var**2)/2
)

show("F2(sum,gap)", F2_sum_gap)
show("expected", expected_sum_gap)

check(
    "sum/gap reconstruction",
    F2_sum_gap,
    expected_sum_gap,
)

print("""
  With

      sum = P+Q
      gap = P-Q

  the layer becomes

      F2 = (sum^2 + 2 sum - 3 gap^2)/2.

  Therefore the F2 layer separates into:

      symmetric channel:
          sum = P+Q = S

      antisymmetric channel:
          gap = P-Q

  without introducing S or Delta as primitive inputs.
""")

print()


# ==============================================================================
# [4] F2 AS A QUADRATIC IN P
# ==============================================================================

print("[4] F2 AS A QUADRATIC IN P")
print("-" * 78)

F2_as_P = sp.Poly(
    sp.expand(F2_PQ - sp.Symbol("F")),
    P,
)

# Rebuild explicitly with F as the observed F2 value.
Fobs = sp.symbols("Fobs")

eq_P = clean(
    F2_PQ - Fobs
)

poly_P = sp.Poly(
    sp.expand(eq_P),
    P,
)

aP = poly_P.coeff_monomial(P**2)
bP = poly_P.coeff_monomial(P)
cP = poly_P.coeff_monomial(1)

show("quadratic coefficient a_P", aP)
show("linear coefficient b_P", bP)
show("constant coefficient c_P", cP)

check("a_P=-1", aP, -1)
check("b_P=4Q+1", bP, 4*Q + 1)
check("c_P=-Q^2+Q-F2", cP, -Q**2 + Q - Fobs)

print()


# ==============================================================================
# [5] DISCRIMINANT OF THE P-QUADRATIC
# ==============================================================================

print("[5] DISCRIMINANT OF THE P-QUADRATIC")
print("-" * 78)

disc_P = clean(
    sp.discriminant(
        -P**2 + (4*Q + 1)*P + (-Q**2 + Q - Fobs),
        P,
    )
)

disc_P_expected = clean(
    12*Q**2 + 12*Q + 1 + 4*Fobs
)

show("disc_P", disc_P)
show("expected", disc_P_expected)

check(
    "discriminant_P",
    disc_P,
    disc_P_expected,
)

print()


# ==============================================================================
# [6] SUBSTITUTE THE ACTUAL F2 EXPRESSION
# ==============================================================================

print("[6] DISCRIMINANT AFTER SUBSTITUTING THE ACTUAL F2")
print("-" * 78)

disc_P_actual = clean(
    disc_P.subs(
        Fobs,
        F2_PQ,
    )
)

disc_P_actual_factored = clean(
    disc_P_actual
)

show(
    "disc_P(F2)",
    disc_P_actual_factored,
)

check(
    "disc_P(F2)=12Q^2+12Q+1+4F2",
    disc_P_actual,
    12*Q**2 + 12*Q + 1 + 4*F2_PQ,
)

print("""
  This is intentionally not simplified away.

  The quantity

      D_P = 12 Q^2 + 12 Q + 1 + 4 F2

  is the exact discriminant of the F2 equation when solving
  for P with Q held fixed.

  The next question is whether this discriminant has a useful
  representation in the original conic coordinates.
""")

print()


# ==============================================================================
# [7] SUBSTITUTE F2 = 2N+S-DELTA
# ==============================================================================

print("[7] DISCRIMINANT IN N,S,DELTA")
print("-" * 78)

disc_P_NSD = clean(
    disc_P.subs(
        Fobs,
        2*N + S - Delta,
    )
)

disc_P_NSD_root = clean(
    disc_P_NSD.subs(Q, Qxy)
)

show("disc_P(N,S,Delta,Q)", disc_P_NSD)
show("disc_P after Q=B/r", disc_P_NSD_root)

# Expected direct P,Q simplification remains the strongest form.
check(
    "discriminant N,S,Delta consistency",
    disc_P_NSD,
    12*Q**2 + 12*Q + 1 + 4*(2*N + S - Delta),
)

print()


# ==============================================================================
# [8] SYMMETRIC COUNTERPART: QUADRATIC IN Q
# ==============================================================================

print("[8] SYMMETRIC COUNTERPART: QUADRATIC IN Q")
print("-" * 78)

eq_Q = clean(
    F2_PQ - Fobs
)

poly_Q = sp.Poly(
    sp.expand(eq_Q),
    Q,
)

aQ = poly_Q.coeff_monomial(Q**2)
bQ = poly_Q.coeff_monomial(Q)
cQ = poly_Q.coeff_monomial(1)

disc_Q = clean(
    sp.discriminant(
        -Q**2 + (4*P + 1)*Q + (-P**2 + P - Fobs),
        Q,
    )
)

expected_disc_Q = clean(
    12*P**2 + 12*P + 1 + 4*Fobs
)

show("a_Q", aQ)
show("b_Q", bQ)
show("c_Q", cQ)
show("disc_Q", disc_Q)

check("a_Q=-1", aQ, -1)
check("b_Q=4P+1", bQ, 4*P + 1)
check("c_Q=-P^2+P-F2", cQ, -P**2 + P - Fobs)
check("disc_Q", disc_Q, expected_disc_Q)

print()


# ==============================================================================
# [9] P <-> Q SYMMETRY
# ==============================================================================

print("[9] P <-> Q SYMMETRY")
print("-" * 78)

check(
    "F2 symmetric under P<->Q",
    F2_PQ,
    F2_PQ.subs({P: Q, Q: P}, simultaneous=True),
)

check(
    "disc_P transforms to disc_Q",
    disc_P.subs({
        Q: P,
    }),
    disc_Q,
)

print()


# ==============================================================================
# [10] DIRECT HISTORICAL x,y CROSSWALK
# ==============================================================================

print("[10] DIRECT HISTORICAL x,y CROSSWALK")
print("-" * 78)

A0 = clean(
    2*y - 2*x + 3
)

B0 = clean(
    2*y + 2*x - 3
)

Pxy_hist = clean(A0 / r)
Qxy_hist = clean(B0 / r)

Nxy = clean(
    Pxy_hist * Qxy_hist
)

Sxy = clean(
    Pxy_hist + Qxy_hist
)

Delta_xy = clean(
    (Pxy_hist - Qxy_hist)**2
)

F2xy = clean(
    F2_PQ.subs({
        P: Pxy_hist,
        Q: Qxy_hist,
    })
)

show("A0", A0)
show("B0", B0)
show("P(x,y)", Pxy_hist)
show("Q(x,y)", Qxy_hist)
show("S(x,y)", Sxy)
show("N(x,y)", Nxy)
show("Delta(x,y)", Delta_xy)
show("F2(x,y)", F2xy)

check(
    "S(x,y)=4y/r",
    Sxy,
    4*y/r,
)

check(
    "N(x,y)=A0B0/r^2",
    Nxy,
    A0*B0/r**2,
)

check(
    "Delta(x,y)=(A0-B0)^2/r^2",
    Delta_xy,
    (A0 - B0)**2/r**2,
)

print()


# ==============================================================================
# [11] F2 AS A QUADRATIC IN y
# ==============================================================================

print("[11] F2 AS A QUADRATIC IN y")
print("-" * 78)

poly_y = sp.Poly(
    sp.expand(F2xy),
    y,
)

ay = clean(poly_y.coeff_monomial(y**2))
by = clean(poly_y.coeff_monomial(y))
cy = clean(poly_y.coeff_monomial(1))

disc_y = clean(
    sp.discriminant(
        sp.expand(F2xy),
        y,
    )
)

show("coefficient y^2", ay)
show("coefficient y", by)
show("constant", cy)
show("discriminant in y", disc_y)

print("""
  This is the key historical-coordinate test:

      S = 4y/r

  is already explicit.

  Therefore F2 should be interpreted as a quadratic observable
  in the original 2020 moving coordinate y, rather than as an
  N-only object.
""")

print()


# ==============================================================================
# [12] F2 AS A QUADRATIC IN x
# ==============================================================================

print("[12] F2 AS A QUADRATIC IN x")
print("-" * 78)

poly_x = sp.Poly(
    sp.expand(F2xy),
    x,
)

ax = clean(poly_x.coeff_monomial(x**2))
bx = clean(poly_x.coeff_monomial(x))
cx = clean(poly_x.coeff_monomial(1))

disc_x = clean(
    sp.discriminant(
        sp.expand(F2xy),
        x,
    )
)

show("coefficient x^2", ax)
show("coefficient x", bx)
show("constant", cx)
show("discriminant in x", disc_x)

print()


# ==============================================================================
# [13] ROOT-ORBIT TRANSLATION
# ==============================================================================

print("[13] ROOT-ORBIT TRANSLATION")
print("-" * 78)

hh = sp.symbols("hh")

A_h = clean(A0 + hh*r)
B_h = clean(B0 + hh*r)

P_h = clean(A_h / r)
Q_h = clean(B_h / r)

N_h = clean(P_h * Q_h)
S_h = clean(P_h + Q_h)
Delta_h = clean((P_h - Q_h)**2)

F2_h = clean(
    F2_PQ.subs({
        P: P_h,
        Q: Q_h,
    })
)

check(
    "N_h=N+hh*S+hh^2",
    N_h,
    Nxy + hh*Sxy + hh**2,
)

check(
    "S_h=S+2hh",
    S_h,
    Sxy + 2*hh,
)

check(
    "Delta_h=Delta",
    Delta_h,
    Delta_xy,
)

check(
    "F2_h-F2=2hh(S+hh+1)",
    F2_h - F2xy,
    2*hh*(Sxy + hh + 1),
)

print()


# ==============================================================================
# [14] F2 ORBIT DISCRIMINANT
# ==============================================================================

print("[14] F2 ORBIT DISCRIMINANT")
print("-" * 78)

F2_orbit_disc = clean(
    sp.discriminant(
        sp.expand(F2_h),
        hh,
    )
)

show("disc_h(T_h(F2))", F2_orbit_disc)

check(
    "disc_h(F2)=12Delta+4",
    F2_orbit_disc,
    12*Delta_xy + 4,
)

print()


# ==============================================================================
# [15] SPECTRAL QUADRATIC DIRECTLY FROM x,y
# ==============================================================================

print("[15] SPECTRAL QUADRATIC DIRECTLY FROM x,y")
print("-" * 78)

Qxy_poly = clean(
    z**2 - Sxy*z + Nxy
)

Qxy_factor = clean(
    (z - Pxy_hist) * (z - Qxy_hist)
)

show("Q_xy(z)", Qxy_poly)
show("factorized Q_xy(z)", Qxy_factor)

check(
    "Q_xy=(z-A0/r)(z-B0/r)",
    Qxy_poly,
    Qxy_factor,
)

check(
    "disc(Q_xy)=Delta",
    sp.discriminant(Qxy_poly, z),
    Delta_xy,
)

print()


# ==============================================================================
# [16] SPECTRAL TRANSLATION
# ==============================================================================

print("[16] SPECTRAL TRANSLATION")
print("-" * 78)

Qxy_shift = clean(
    Qxy_poly.subs(z, z - 1)
)

Qxy_shift_expected = clean(
    (z - Pxy_hist - 1)
    * (z - Qxy_hist - 1)
)

check(
    "Q_xy(z-1) translated roots",
    Qxy_shift,
    Qxy_shift_expected,
)

print()


# ==============================================================================
# [17] ROOT ORBIT REFLECTION
# ==============================================================================

print("[17] ROOT ORBIT REFLECTION")
print("-" * 78)

L_h = clean(
    (A0 + hh*r) * (B0 + hh*r) / r**2
)

check(
    "L_h=N+hh*S+hh^2",
    L_h,
    Nxy + hh*Sxy + hh**2,
)

check(
    "L_-z=Q_xy(z)",
    L_h.subs(hh, -z),
    Qxy_poly,
)

check(
    "L(h+2)-2L(h+1)+L(h)=2",
    L_h.subs(hh, hh+2)
    - 2*L_h.subs(hh, hh+1)
    + L_h,
    2,
)

print()


# ==============================================================================
# [18] DIRECT ROOT-COORDINATE SUMMARY
# ==============================================================================

print("[18] DIRECT ROOT-COORDINATE SUMMARY")
print("-" * 78)

print(
r"""
  The 2020 coordinates give the complete moving/invariant split:

      A0 = 2y - 2x + 3
      B0 = 2y + 2x - 3

      A0+B0 = 4y
      B0-A0 = 4x-6
      A0B0   = 4v-1

  After normalization:

      P = A0/r
      Q = B0/r

      S     = P+Q
      N     = PQ
      Delta = (P-Q)^2

  The known homogeneous layer is:

      F2 = 4PQ-P^2-Q^2+P+Q

         = 2N+S-Delta.

  The root translation is:

      A0 -> A0 + h r
      B0 -> B0 + h r

  giving:

      L_h = N+hS+h^2.

  The spectral polynomial is the reflection:

      Q(z)=L(-z).

  This experiment intentionally does NOT assume anything about
  F4, F5, or F6.
"""
)

print()


# ==============================================================================
# [19] FINAL AUDIT
# ==============================================================================

print("[19] FINAL AUDIT")
print("-" * 78)

print(f"  built-in symbolic failures = {failures}")
print("  higher homogeneous layers = NOT ASSUMED")
print("  F2/F3 root-coordinate structure = TESTED")

print()
print("=" * 78)
print("EXPERIMENT 564 FINISHED")
print("=" * 78)
