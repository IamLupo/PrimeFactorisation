#!/usr/bin/env python3

import sympy as sp


print("=" * 78)
print("EXPERIMENT 517 START")
print("=" * 78)
print("MINIMAL KAPPA OBSERVABLE -> RECURRENCE -> S, DELTA")
print("=" * 78)


# ============================================================================
# Helpers
# ============================================================================

def fact(x):
    return sp.factor(sp.expand(x))


def simp(x):
    return sp.simplify(sp.expand(x))


def zero(x):
    return simp(x) == 0


def cert(name, expr):
    d = sp.factor(sp.simplify(expr))
    ok = d == 0
    print(f"  {name}")
    print(f"    difference = {d}")
    print(f"    PASS = {ok}")
    return ok


# ============================================================================
# KAPPA sequence
# ============================================================================

P, Q = sp.symbols("P Q")
N, S = sp.symbols("N S")


def F(n, P_, Q_):
    return sp.expand(
        P_ * (Q_ + 1) ** n
        + Q_ * (P_ + 1) ** n
        - (Q_ + 1) * P_ ** n
        - (P_ + 1) * Q_ ** n
    )


def Fn_int(n, p, q):
    return (
        p * (q + 1) ** n
        + q * (p + 1) ** n
        - (q + 1) * p ** n
        - (p + 1) * q ** n
    )


# ============================================================================
# Prime test set
# ============================================================================

PRIME_PAIRS = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),
    (3, 5),
    (5, 13),
    (13, 17),
    (3, 7),
    (7, 11),
]


DEGENERATE_PAIRS = [
    (2, 3),
    (3, 4),
    (4, 5),
    (5, 6),
    (3, 6),
    (4, 8),
]


# ============================================================================
# [1] SPECTRAL COLLISION AUDIT
# ============================================================================

print()
print("[1] SPECTRAL BASE COLLISION AUDIT")
print("-" * 78)

def bases(p, q):
    return [p, q, p + 1, q + 1]


def spectral_rank(p, q):
    vals = bases(p, q)
    return len(set(vals))


for p, q in PRIME_PAIRS + DEGENERATE_PAIRS:
    vals = bases(p, q)
    rank = len(set(vals))
    print(
        f"  ({p},{q})"
        f" bases={vals}"
        f" distinct={rank}"
        f" generic_order4={rank == 4}"
    )


# ============================================================================
# [2] SYMBOLIC RECURRENCE
# ============================================================================

print()
print("[2] SYMBOLIC ORDER-4 RECURRENCE")
print("-" * 78)

z = sp.symbols("z")

chi = sp.expand(
    (z - P)
    * (z - Q)
    * (z - P - 1)
    * (z - Q - 1)
)

poly = sp.Poly(chi, z)

a1 = sp.expand(-poly.all_coeffs()[1])
a2 = sp.expand(poly.all_coeffs()[2])
a3 = sp.expand(-poly.all_coeffs()[3])
a4 = sp.expand(poly.all_coeffs()[4])

print(f"  a1 = {fact(a1)}")
print(f"  a2 = {fact(a2)}")
print(f"  a3 = {fact(a3)}")
print(f"  a4 = {fact(a4)}")


# ============================================================================
# [3] S AND DELTA EXTRACTION FROM a1,a2
# ============================================================================

print()
print("[3] RECURRENCE -> SYMMETRIC COORDINATES")
print("-" * 78)

a1_target = 2 * (P + Q + 1)
a2_target = P**2 + 4*P*Q + 3*P + Q**2 + 3*Q + 1

S_from_a1 = sp.expand(a1 / 2 - 1)

Delta_from_a1a2 = sp.expand(
    sp.Rational(3, 4) * a1**2
    - 2 * a2
    - 1
)

delta_target = sp.expand((P - Q)**2)

cert(
    "S = a1/2 - 1",
    S_from_a1 - (P + Q),
)

cert(
    "Delta = 3a1^2/4 - 2a2 - 1",
    Delta_from_a1a2 - delta_target,
)


# ============================================================================
# [4] HANKEL MATRIX FOR THE KAPPA SEQUENCE
# ============================================================================
#
# For a sequence satisfying
#
#   F[n+4] = a1 F[n+3] - a2 F[n+2] + a3 F[n+1] - a4 F[n],
#
# define
#
#   H4(n) = det(F[n+i+j])_{0<=i,j<4}.
#
# The determinant vanishes exactly when the spectral rank is < 4.
#
# We test that relation explicitly.
# ============================================================================

print()
print("[4] HANKEL RANK / SPECTRAL DEGENERACY")
print("-" * 78)


def hankel(values, start, size):
    return sp.Matrix([
        [
            values[start + i + j]
            for j in range(size)
        ]
        for i in range(size)
    ])


hankel_failures = 0

for p, q in PRIME_PAIRS + DEGENERATE_PAIRS:

    vals = [
        sp.Integer(Fn_int(n, p, q))
        for n in range(12)
    ]

    H4 = sp.factor(hankel(vals, 0, 4).det())
    H5 = sp.factor(
        sp.Matrix([
            [vals[i + j] for j in range(5)]
            for i in range(5)
        ]).det()
    )

    distinct = len(set(bases(p, q)))

    print(
        f"  ({p},{q})"
        f" distinct_bases={distinct}"
        f" H4={H4}"
        f" H5={H5}"
    )

    expected_rank4 = distinct == 4

    if expected_rank4 and H4 == 0:
        hankel_failures += 1

    if distinct < 5 and H5 != 0:
        # Any four-base sequence must have H5 = 0.
        hankel_failures += 1

print()
print(f"  HANKEL FAILURES = {hankel_failures}")


# ============================================================================
# [5] MINIMAL TERM COUNT FOR FOURTH-ORDER RECURRENCE
# ============================================================================
#
# Four unknown coefficients require four recurrence equations.
#
# Each equation needs F[n],...,F[n+4].
#
# Therefore F_0...F_7 is theoretically the minimal generic window.
#
# We test:
#
#   F_0...F_7
#   F_1...F_8
#   F_2...F_9
#
# and verify that all produce the same recurrence.
# ============================================================================

print()
print("[5] MINIMAL OBSERVABLE WINDOW")
print("-" * 78)


def recover_order4(values, offset=0):
    rows = []
    rhs = []

    for n in range(offset, offset + 4):
        rows.append([
            -values[n + 3],
             values[n + 2],
            -values[n + 1],
             values[n],
        ])
        rhs.append(-values[n + 4])

    A = sp.Matrix(rows)
    b = sp.Matrix(rhs)

    if A.det() == 0:
        return None

    sol = A.LUsolve(b)

    return tuple(sp.simplify(x) for x in sol)


minimal_failures = 0

for p, q in PRIME_PAIRS:

    vals = [
        sp.Integer(Fn_int(n, p, q))
        for n in range(10)
    ]

    rec0 = recover_order4(vals, 0)
    rec1 = recover_order4(vals, 1)
    rec2 = recover_order4(vals, 2)

    consistent = (
        rec0 is not None
        and rec1 is not None
        and rec2 is not None
        and rec0 == rec1 == rec2
    )

    print(
        f"  ({p},{q})"
        f" window0={rec0}"
        f" window1={rec1}"
        f" window2={rec2}"
        f" stable={consistent}"
    )

    if not consistent:
        minimal_failures += 1

print()
print(f"  MINIMAL WINDOW FAILURES = {minimal_failures}")


# ============================================================================
# [6] DIRECT RECURRENCE COEFFICIENT FORMULAS FROM HANKEL MINORS
# ============================================================================
#
# Cramer's rule expresses the recurrence coefficients directly as
# rational functions of sequence values.
#
# This section constructs those formulas symbolically.
# ============================================================================

print()
print("[6] RECURRENCE COEFFICIENTS AS HANKEL-MINOR RATIOS")
print("-" * 78)

f = sp.symbols("f0:8")

# Four equations:
#
# -a1 f[n+3] + a2 f[n+2] - a3 f[n+1] + a4 f[n] = -f[n+4]

Arec = sp.Matrix([
    [-f[3], f[2], -f[1], f[0]],
    [-f[4], f[3], -f[2], f[1]],
    [-f[5], f[4], -f[3], f[2]],
    [-f[6], f[5], -f[4], f[3]],
])

brec = sp.Matrix([
    -f[4],
    -f[5],
    -f[6],
    -f[7],
])

detA = sp.factor(Arec.det())

print(f"  recurrence matrix determinant =")
print(f"    {detA}")

if detA != 0:
    rec_symbols = sp.Matrix(Arec).inv() * brec

    R1 = sp.factor(rec_symbols[0])
    R2 = sp.factor(rec_symbols[1])
    R3 = sp.factor(rec_symbols[2])
    R4 = sp.factor(rec_symbols[3])

    print()
    print(f"  a1(F0..F7) = {R1}")
    print()
    print(f"  a2(F0..F7) = {R2}")
    print()
    print(f"  a3(F0..F7) = {R3}")
    print()
    print(f"  a4(F0..F7) = {R4}")
else:
    R1 = R2 = R3 = R4 = None


# ============================================================================
# [7] FORMAL OBSERVABLE GAP FUNCTION
# ============================================================================

print()
print("[7] FORMAL OBSERVABLE GAP FUNCTION")
print("-" * 78)

if R1 is not None:

    Delta_observable = sp.factor(
        sp.Rational(3, 4) * R1**2
        - 2 * R2
        - 1
    )

    S_observable = sp.factor(
        R1 / 2 - 1
    )

    print("  S(F0..F7) =")
    print(f"    {S_observable}")

    print()
    print("  Delta(F0..F7) =")
    print(f"    {Delta_observable}")

    print()
    print("  This is an explicit rational function of the")
    print("  eight observable KAPPA terms.")


# ============================================================================
# [8] NUMERICAL DIRECT F-ONLY EXTRACTION
# ============================================================================
#
# IMPORTANT:
# The calculation below intentionally knows only F_0...F_7.
# p and q are used solely to manufacture the test sequence.
# The extraction itself does not receive them.
# ============================================================================

print()
print("[8] F-ONLY RECOVERY AUDIT")
print("-" * 78)

f_only_failures = 0

for p, q in PRIME_PAIRS:

    observed = [
        sp.Integer(Fn_int(n, p, q))
        for n in range(8)
    ]

    rec = recover_order4(observed, 0)

    if rec is None:
        print(f"  ({p},{q}) recurrence unavailable")
        f_only_failures += 1
        continue

    rr1, rr2, rr3, rr4 = rec

    S_rec = sp.simplify(rr1 / 2 - 1)
    Delta_rec = sp.simplify(
        sp.Rational(3, 4) * rr1**2
        - 2 * rr2
        - 1
    )

    N_value = p * q

    roots = sp.solve(
        sp.Symbol("z")**2 - S_rec*sp.Symbol("z") + N_value,
        sp.Symbol("z"),
    )

    expected_roots = {sp.Integer(p), sp.Integer(q)}

    # solve() may return symbolic ordering.
    roots_set = set(sp.simplify(r) for r in roots)

    okS = S_rec == p + q
    okD = Delta_rec == (p - q) ** 2
    okRoots = roots_set == expected_roots

    ok = okS and okD and okRoots

    print(
        f"  ({p},{q})"
        f" S={S_rec}"
        f" Delta={Delta_rec}"
        f" roots={roots_set}"
        f" PASS={ok}"
    )

    if not ok:
        f_only_failures += 1

print()
print(f"  F-ONLY FAILURES = {f_only_failures}")


# ============================================================================
# [9] DEGENERATE SPECTRAL CASES
# ============================================================================

print()
print("[9] DEGENERATE SPECTRAL CASES")
print("-" * 78)

degenerate_failures = 0

for p, q in DEGENERATE_PAIRS:

    values = [
        sp.Integer(Fn_int(n, p, q))
        for n in range(12)
    ]

    distinct = len(set(bases(p, q)))

    rec = recover_order4(values, 0)

    # Expected:
    #   unique order-4 recovery only when all four bases are distinct.
    #
    # When collision occurs, failure is acceptable and informative.

    expected_generic = distinct == 4

    acceptable = (
        rec is not None
        if expected_generic
        else True
    )

    print(
        f"  ({p},{q})"
        f" bases={bases(p,q)}"
        f" distinct={distinct}"
        f" order4_recovery={rec}"
        f" ACCEPTED={acceptable}"
    )

    if not acceptable:
        degenerate_failures += 1

print()
print(f"  DEGENERACY FAILURES = {degenerate_failures}")


# ============================================================================
# [10] PRIME DOMAIN NONDEGENERACY THEOREM
# ============================================================================

print()
print("[10] PRIME-DOMAIN NONDEGENERACY")
print("-" * 78)

print("""
For distinct odd primes p,q:

    P != Q
    P != Q+1
    Q != P+1

because two distinct odd primes cannot differ by 1.

Therefore

    {P,Q,P+1,Q+1}

contains four distinct values.

Hence the order-4 recurrence is generically identifiable
on the intended distinct-odd-prime domain.

The failures seen in Experiment 516 were caused by controls
such as (3,4), (4,5), and (5,6), which are outside that domain.
""")


# ============================================================================
# [11] INFORMATION-BOUNDARY CHECK
# ============================================================================

print()
print("[11] INFORMATION-BOUNDARY CHECK")
print("-" * 78)

print("""
The downstream mechanism is now explicit:

    F0,...,F7
       |
       v
    recurrence coefficients
       |
       +----> a1
       |       |
       |       v
       |      S = a1/2 - 1
       |
       +----> a2
               |
               v
        Delta = 3a1^2/4 - 2a2 - 1
               |
               v
          z^2 - S z + N

This is not yet an N-only factorization algorithm.

The unresolved interface remains:

    N
      -> independently generated F0,...,F7
      -> recurrence
      -> S, Delta
      -> p,q.

Therefore the next useful question is extremely specific:

    Can the original homogeneous-layer construction generate
    the eight F-values, or enough equivalent invariants to
    recover the same recurrence matrix, without already using
    p or q?
""")


# ============================================================================
# [12] FINAL AUDIT
# ============================================================================

print()
print("[12] FINAL EXACT AUDIT")
print("-" * 78)

checks = [
    ("symbolic Delta identity", zero(
        sp.Rational(3, 4) * a1**2
        - 2*a2
        - 1
        - (P-Q)**2
    )),
    ("symbolic recurrence", all(
        zero(
            F(n+4, P, Q)
            - a1*F(n+3, P, Q)
            + a2*F(n+2, P, Q)
            - a3*F(n+1, P, Q)
            + a4*F(n, P, Q)
        )
        for n in range(8)
    )),
    ("Hankel audit", hankel_failures == 0),
    ("minimal windows", minimal_failures == 0),
    ("F-only recovery", f_only_failures == 0),
    ("degeneracy handling", degenerate_failures == 0),
]

for name, ok in checks:
    print(f"  {name:<28} = {ok}")

overall = all(ok for _, ok in checks)

print()
print(f"  OVERALL EXACT AUDIT = {overall}")


# ============================================================================
# [13] CONCLUSION
# ============================================================================

print()
print("[13] RESEARCH CONCLUSION")
print("-" * 78)

print(
r"""
Experiment 517 establishes the minimal downstream observable
interface for the generic four-base KAPPA sequence.

For distinct odd primes:

    F0,...,F7

are sufficient to reconstruct the order-4 recurrence

    F[n+4]
      = a1 F[n+3]
      - a2 F[n+2]
      + a3 F[n+1]
      - a4 F[n].

From only a1 and a2:

    S
      = a1/2 - 1,

and

    Delta
      = 3a1^2/4 - 2a2 - 1.

Therefore:

    F0,...,F7
       -> a1,a2
       -> S,Delta
       -> z^2-Sz+N
       -> p,q.

The spectral-degeneracy controls show why non-prime examples
can fail: repeated bases reduce the true recurrence order.

The next research step should therefore move upstream.

Do NOT generate more identities involving S or Delta.

Instead inspect the existing N-only homogeneous-layer objects
and ask whether eight quantities equivalent to

    F0,...,F7

or four independent recurrence equations

can be obtained without constructing p,q.

That is now the narrowest remaining bridge.
"""
)

print()
print("=" * 78)
print("EXPERIMENT 517 FINISHED")
print("=" * 78)
