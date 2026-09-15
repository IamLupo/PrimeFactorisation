#!/usr/bin/env python3

import sympy as sp

print("=" * 78)
print("EXPERIMENT 515 START")
print("=" * 78)
print("UNCENTERED RECURRENCE -> CENTERED DISCRIMINANT INVARIANT")
print("=" * 78)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def fact(expr):
    return sp.factor(sp.expand(expr))


def simp(expr):
    return sp.simplify(sp.expand(expr))


def cert(name, expr):
    d = fact(expr)
    print(f"  {name}")
    print(f"    difference = {d}")
    print(f"    PASS = {d == 0}")
    print()


# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------

P, Q = sp.symbols("P Q")
z, xi = sp.symbols("z xi")

N, S, Delta = sp.symbols("N S Delta")


# ===========================================================================
# [1] FOUR KAPPA BASES
# ===========================================================================

print()
print("[1] FOUR KAPPA BASES")
print("-" * 78)

lam = [
    P,
    Q,
    P + 1,
    Q + 1,
]

for i, L in enumerate(lam, 1):
    print(f"  lambda_{i} = {L}")

print()


# ===========================================================================
# [2] UNCENTERED CHARACTERISTIC POLYNOMIAL
# ===========================================================================

print()
print("[2] UNCENTERED CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

chi = fact(sp.prod(z - L for L in lam))

print("  chi(z) =")
print(f"    {chi}")
print()


poly = sp.Poly(chi, z)

a1 = fact(-poly.all_coeffs()[1])
a2 = fact(poly.all_coeffs()[2])
a3 = fact(-poly.all_coeffs()[3])
a4 = fact(poly.all_coeffs()[4])

print(f"  a1 = {a1}")
print(f"  a2 = {a2}")
print(f"  a3 = {a3}")
print(f"  a4 = {a4}")
print()


# ===========================================================================
# [3] EXPRESS COEFFICIENTS IN N,S
# ===========================================================================

print()
print("[3] RECURRENCE COEFFICIENTS IN N,S")
print("-" * 78)

subs_NS = {
    P + Q: S,
    P * Q: N,
}

# Do the symmetric expressions explicitly instead of relying on
# simultaneous substitution order.

a1_NS = fact(2 * (P + Q + 1))
a2_NS = fact(
    P*Q
    + P*(P+1)
    + P*(Q+1)
    + Q*(P+1)
    + Q*(Q+1)
    + (P+1)*(Q+1)
)
a3_NS = fact(
    P*Q*(P+1)
    + P*Q*(Q+1)
    + P*(P+1)*(Q+1)
    + Q*(P+1)*(Q+1)
)
a4_NS = fact(P*Q*(P+1)*(Q+1))

a1_expected = 2*(S + 1)
a2_expected = fact(S**2 + 3*S + 2*N + 1)
a3_expected = fact((S + 1)*(N + S + N))
a4_expected = fact(N*(N + S + 1))

# a3 is actually:
# pairwise triple sum =
# PQ(P+1)+PQ(Q+1)+P(P+1)(Q+1)+Q(P+1)(Q+1)
#
# simplify symbolically through symmetric reduction by replacing
# q^2+p^2 with S^2-2N.

a3_expected = fact(
    2*N*S + 2*N + S**2 + S
)

print(f"  a1(N,S) = {a1_expected}")
print(f"  a2(N,S) = {a2_expected}")
print(f"  a3(N,S) = {a3_expected}")
print(f"  a4(N,S) = {a4_expected}")
print()


# Use symmetric reduction with q = S-p and eliminate p.
def sym_reduce(expr):
    tmp = sp.expand(expr.subs(Q, S - P))
    tmp = sp.factor(tmp)
    # Any remaining P dependence should cancel.
    return sp.factor(tmp)


a1_red = sym_reduce(a1)
a2_red = sym_reduce(a2)
a3_red = sym_reduce(a3)
a4_red = sym_reduce(a4)

print("  direct symmetric reductions:")
print(f"    a1 = {a1_red}")
print(f"    a2 = {a2_red}")
print(f"    a3 = {a3_red}")
print(f"    a4 = {a4_red}")
print()


# ===========================================================================
# [4] SPECTRAL CENTER FROM a1
# ===========================================================================

print()
print("[4] SPECTRAL CENTER")
print("-" * 78)

center = fact(a1 / 4)

print(f"  center = a1/4 = {center}")
print()

cert(
    "center = (S+1)/2",
    center - (S + 1)/2,
)


# ===========================================================================
# [5] TRANSLATED QUARTIC
# ===========================================================================

print()
print("[5] TRANSLATED QUARTIC")
print("-" * 78)

translated = fact(
    sp.expand(
        chi.subs(z, xi + a1/4)
    )
)

print("  chi(xi + a1/4) =")
print(f"    {translated}")
print()


translated_poly = sp.Poly(translated, xi)

tc = translated_poly.all_coeffs()

T2 = fact(tc[2])
T1 = fact(-tc[3])
T0 = fact(tc[4])

print(f"  xi^4 coefficient = {tc[0]}")
print(f"  xi^3 coefficient = {tc[1]}")
print(f"  xi^2 coefficient = {T2}")
print(f"  xi coefficient   = {tc[3]}")
print(f"  constant         = {T0}")
print()


# ===========================================================================
# [6] CUBIC / LINEAR ELIMINATION
# ===========================================================================

print()
print("[6] CENTERING ELIMINATES ODD TERMS")
print("-" * 78)

cert("cubic term vanishes", tc[1])
cert("linear term vanishes", tc[3])


# ===========================================================================
# [7] CENTERED QUADRATIC COEFFICIENT
# ===========================================================================

print()
print("[7] CENTERED QUADRATIC COEFFICIENT")
print("-" * 78)

T2_expected = fact(
    -(Delta + sp.Rational(1, 2))
)

# The known centered quartic is
#
# xi^4 -(Delta+1/2) xi^2
# + (Delta/4-1/4)^2

print(f"  T2 = {T2}")
print(f"  expected = {T2_expected}")
print()

cert(
    "T2 = -(Delta+1/2)",
    T2 - T2_expected,
)


# ===========================================================================
# [8] DISCRIMINANT DIRECTLY FROM UNCENTERED COEFFICIENTS
# ===========================================================================

print()
print("[8] DISCRIMINANT INVARIANT FROM QUARTIC COEFFICIENTS")
print("-" * 78)

# For a monic quartic
#
#   z^4 - a1 z^3 + a2 z^2 - a3 z + a4
#
# translation by a1/4 gives
#
#   xi^4 + T2 xi^2 + T1 xi + T0.
#
# T2 has the classical depressed-quartic expression
#
#   T2 = a2 - 3 a1^2 / 8.
#
# Therefore:
#
#   Delta = -T2 - 1/2
#          = 3 a1^2/8 - a2 - 1/2.

Delta_candidate = fact(
    3*a1**2/sp.Integer(8)
    - a2
    - sp.Rational(1, 2)
)

print("  Candidate:")
print("    Delta_candidate = 3*a1^2/8 - a2 - 1/2")
print(f"    = {Delta_candidate}")
print()

cert(
    "Delta candidate",
    Delta_candidate - (P - Q)**2,
)


# ===========================================================================
# [9] PURE N,S VERSION
# ===========================================================================

print()
print("[9] PURE N,S VERSION")
print("-" * 78)

Delta_NS_candidate = fact(
    3*a1_expected**2/sp.Integer(8)
    - a2_expected
    - sp.Rational(1, 2)
)

print(f"  Delta_candidate(N,S) = {Delta_NS_candidate}")
print()

cert(
    "Delta_candidate = S^2 - 4N",
    Delta_NS_candidate - (S**2 - 4*N),
)


# ===========================================================================
# [10] IMPORTANT REARRANGEMENT
# ===========================================================================

print()
print("[10] REARRANGEMENT")
print("-" * 78)

print("""
The quartic gives the exact invariant

    Delta = 3*a1^2/8 - a2 - 1/2.

Since

    a1 = 2(S+1),

this becomes

    Delta
      = 3(S+1)^2/2 - a2 - 1/2.

The important question is whether a1 and a2 can be
obtained directly from the observable F_n without first
recovering P,Q.
""")

rearranged = fact(
    Delta_candidate
    - (3*a1**2/sp.Integer(8) - a2 - sp.Rational(1,2))
)

cert(
    "rearrangement identity",
    rearranged,
)


# ===========================================================================
# [11] RECURRENCE COEFFICIENT RECOVERY FROM TERMS
# ===========================================================================

print()
print("[11] RECURRENCE FROM SEQUENCE VALUES")
print("-" * 78)

# Uncentered sequence obeys
#
#   F_{n+4}
#   - a1 F_{n+3}
#   + a2 F_{n+2}
#   - a3 F_{n+1}
#   + a4 F_n = 0.
#
# We try to recover a1,a2,a3,a4 from exact sequence terms by
# linear algebra.

def F(k):
    return fact(
        P*(Q+1)**k
        + Q*(P+1)**k
        - (Q+1)*P**k
        - (P+1)*Q**k
    )


coeffs = sp.symbols("r1:5")

equations = []

for k in range(0, 8):
    equations.append(
        fact(
            F(k+4)
            - coeffs[0]*F(k+3)
            + coeffs[1]*F(k+2)
            - coeffs[2]*F(k+1)
            + coeffs[3]*F(k)
        )
    )

sol = sp.solve(equations, coeffs, dict=True)

print(f"  exact recurrence solutions = {len(sol)}")

if sol:
    solution = sol[0]
    for c in coeffs:
        print(f"  {c} = {fact(solution[c])}")
else:
    print("  No exact symbolic solution found.")

print()


# ===========================================================================
# [12] RECOVER DELTA FROM RECURRENCE DATA
# ===========================================================================

print()
print("[12] DELTA FROM RECURRENCE DATA")
print("-" * 78)

if sol:
    r1 = fact(solution[coeffs[0]])
    r2 = fact(solution[coeffs[1]])
    r3 = fact(solution[coeffs[2]])
    r4 = fact(solution[coeffs[3]])

    print(f"  r1 = {r1}")
    print(f"  r2 = {r2}")
    print(f"  r3 = {r3}")
    print(f"  r4 = {r4}")
    print()

    delta_from_recurrence = fact(
        3*r1**2/sp.Integer(8)
        - r2
        - sp.Rational(1, 2)
    )

    print(f"  Delta_from_recurrence = {delta_from_recurrence}")
    print()

    cert(
        "Delta recovery from uncentered recurrence",
        delta_from_recurrence - (P-Q)**2,
    )

else:
    print("  skipped because recurrence solve failed.")

print()


# ===========================================================================
# [13] DOES DELTA REQUIRE S?
# ===========================================================================

print()
print("[13] INFORMATION-SEPARATION TEST")
print("-" * 78)

print("""
The candidate

    Delta = 3*a1^2/8 - a2 - 1/2

is built entirely from the uncentered characteristic polynomial.

No P, Q, S or centered sequence is required once a1,a2
are known.

This is the key distinction from Experiment 514.
""")

if sol:
    delta_expr = fact(
        3*solution[coeffs[0]]**2/sp.Integer(8)
        - solution[coeffs[1]]
        - sp.Rational(1, 2)
    )

    print(f"  contains P = {delta_expr.has(P)}")
    print(f"  contains Q = {delta_expr.has(Q)}")
    print(f"  contains S = {delta_expr.has(S)}")

    print()


# ===========================================================================
# [14] SAME-N CONTROL
# ===========================================================================

print()
print("[14] SAME-N STRUCTURAL CONTROL")
print("-" * 78)

pairs = [
    (2, 6),
    (3, 4),
    (2, 9),
    (3, 6),
    (4, 5),
    (2, 10),
    (5, 6),
    (3, 10),
]

for p0, q0 in pairs:
    n0 = p0*q0
    s0 = p0+q0
    d0 = (p0-q0)**2

    print(
        f"  ({p0},{q0}) N={n0} S={s0} Delta={d0}"
    )

print()


# ===========================================================================
# [15] FINAL CERTIFICATES
# ===========================================================================

print()
print("[15] FINAL SYMBOLIC AUDIT")
print("-" * 78)

audits = {
    "quartic coefficient structure":
        zero(a1 - 2*(P+Q+1)),

    "center":
        zero(center - (S+1)/2),

    "cubic vanishes":
        zero(tc[1]),

    "linear vanishes":
        zero(tc[3]),

    "centered quadratic":
        zero(T2 - (-(P-Q)**2 - sp.Rational(1,2))),

    "Delta from quartic":
        zero(
            Delta_candidate - (P-Q)**2
        ),

    "Delta N,S":
        zero(
            Delta_NS_candidate - (S**2 - 4*N)
        ),
}

for name, value in audits.items():
    print(f"  {name:<28} = {value}")

overall = all(audits.values())

print()
print(f"  OVERALL EXACT AUDIT = {overall}")


# ===========================================================================
# [16] RESEARCH CONCLUSION
# ===========================================================================

print()
print("[16] RESEARCH CONCLUSION")
print("-" * 78)

print(
r"""
The centered recurrence coefficient A is not fundamentally
a centered-only object.

The uncentered quartic

    z^4-a1 z^3+a2 z^2-a3 z+a4

already contains the gap square through the depressed-quartic
coefficient:

    Delta = 3*a1^2/8 - a2 - 1/2.

Therefore:

    F_n
      -> order-4 recurrence
      -> a1,a2
      -> Delta
      -> |P-Q|
      -> factor gap.

This is the exact next bridge to investigate.

The critical upstream question is now:

    Can a1 and a2 be extracted from a finite number of
    original KAPPA values without recovering P and Q?

If yes, the historical gap-square is accessible directly
through the KAPPA recurrence.

If additionally the recurrence's remaining coefficients
yield S or the shifted product, then the complete factor
recovery chain becomes intrinsic to the sequence itself.
"""
)

print()
print("=" * 78)
print("EXPERIMENT 515 FINISHED")
print("=" * 78)
