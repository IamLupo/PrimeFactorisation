import sympy as sp
from math import isqrt


# ============================================================
# KAPPA CONTROLLED-VARIABLE / FACTOR-RECOVERY ANALYSIS
# ============================================================

print("=" * 80)
print("KAPPA CONTROLLED-VARIABLE / FACTOR-RECOVERY ANALYSIS")
print("=" * 80)


# ============================================================
# SYMBOLS AND BASIC DEFINITIONS
# ============================================================

p, q, r, s = sp.symbols("p q r s")
u, v, A = sp.symbols("u v A")

def F(x):
    return x**2 - x + 1


# ============================================================
# 1. EXACT KAPPA STRUCTURE
# ============================================================

print("\n" + "=" * 80)
print("1. EXACT STRUCTURE")
print("=" * 80)

K3 = sp.expand(1 - F(p) * F(q) * F(r))
K4 = sp.expand(1 - F(p) * F(q) * F(r) * F(s))

print("K3 factorized:")
print(sp.factor(K3))

print("\nK4 factorized:")
print(sp.factor(K4))

transition34 = sp.factor(K4 - K3)

print("\nK4 - K3:")
print(transition34)

expected_transition = sp.expand(
    s * (1 - s) * (1 - K3)
)

residual = sp.factor(
    transition34 - expected_transition
)

print("\nTransition residual:")
print(residual)

assert residual == 0

print("\nPASS: exact 3 -> 4 transition confirmed.")


# ============================================================
# 2. EXACT SYMMETRIC REDUCTION
# ============================================================

print("\n" + "=" * 80)
print("2. SYMMETRIC REDUCTION OF F(p)F(q)")
print("=" * 80)

FpFq = sp.expand(F(p) * F(q))

print("Direct:")
print(FpFq)

# Correct formula:
#
# F(p)F(q)
# = u^2 - uv - u + v^2 - v + 1
#
# where:
#     u = p + q
#     v = pq

FpFq_uv = sp.expand(
    u**2 - u*v - u + v**2 - v + 1
)

print("\nIn terms of u=p+q and v=pq:")
print(FpFq_uv)

verification = sp.expand(
    FpFq -
    FpFq_uv.subs({
        u: p + q,
        v: p * q
    })
)

print("\nSymmetric identity residual:")
print(sp.factor(verification))

assert verification == 0

print("\nPASS: symmetric reduction confirmed.")


# ============================================================
# 3. KAPPA_3 AS A QUADRATIC IN p+q
# ============================================================

print("\n" + "=" * 80)
print("3. KAPPA_3 AND THE p+q RECOVERY EQUATION")
print("=" * 80)

print("""
Let

    A = F(p)F(q)

and

    v = pq = n
    u = p+q.

Then

    A = u^2 - uv - u + v^2 - v + 1.

Therefore

    u^2 -(v+1)u +(v^2-v+1-A) = 0.
""")

recovery_polynomial = sp.expand(
    u**2
    - (v + 1) * u
    + (v**2 - v + 1 - A)
)

print("Recovery polynomial:")
print(recovery_polynomial)

discriminant_u = sp.factor(
    sp.discriminant(recovery_polynomial, u)
)

print("\nDiscriminant for u=p+q:")
print(discriminant_u)

expected_discriminant = sp.expand(
    4*A - 3*(v - 1)**2
)

print("\nSimplified form:")
print(expected_discriminant)

assert sp.expand(
    discriminant_u - expected_discriminant
) == 0

print("\nPASS: recovery discriminant confirmed.")


# ============================================================
# 4. COMPLETE SYMBOLIC FACTOR-RECOVERY CHAIN
# ============================================================

print("\n" + "=" * 80)
print("4. COMPLETE SYMBOLIC FACTOR-RECOVERY CHAIN")
print("=" * 80)

print("""
Controlled r gives:

    1 - K3 = F(p)F(q)F(r)

so

    A = (1-K3)/F(r).

Since n=pq:

    D_u = 4A - 3(n-1)^2

and therefore

    p+q = (n+1 +/- sqrt(D_u))/2.

After recovering u=p+q:

    p,q are roots of

        x^2 - u*x + n = 0.

Thus

    D_factor = u^2 - 4n

and

    p,q = (u +/- sqrt(D_factor))/2.
""")


# ============================================================
# 5. NUMERICAL FACTOR RECOVERY FUNCTION
# ============================================================

def recover_factors(n, kappa3, r):
    """
    Recover p,q from:

        n = p*q
        kappa3 = 1 - F(p)F(q)F(r)

    assuming r is known.
    """

    fr = int(F(r))
    numerator = int(1 - kappa3)

    result = {
        "success": False,
        "n": n,
        "r": r,
        "F_r": fr,
        "A": None,
        "D_u": None,
        "u_candidates": [],
        "factor_candidates": []
    }

    if fr == 0:
        result["reason"] = "F(r)=0"
        return result

    if numerator % fr != 0:
        result["reason"] = "F(r) does not divide 1-kappa3"
        return result

    A_value = numerator // fr
    result["A"] = A_value

    # D_u = 4A - 3(n-1)^2
    D_u = 4*A_value - 3*(n - 1)**2
    result["D_u"] = D_u

    if D_u < 0:
        result["reason"] = "negative discriminant for p+q"
        return result

    sqrt_Du = isqrt(D_u)

    if sqrt_Du**2 != D_u:
        result["reason"] = "D_u is not a perfect square"
        return result

    # u = (n+1 +/- sqrt(D_u))/2
    for sign in (+1, -1):

        numerator_u = n + 1 + sign * sqrt_Du

        if numerator_u % 2 != 0:
            continue

        u_value = numerator_u // 2

        result["u_candidates"].append(u_value)

        # Now solve x^2-u*x+n=0
        D_factor = u_value**2 - 4*n

        if D_factor < 0:
            continue

        sqrt_Dfactor = isqrt(D_factor)

        if sqrt_Dfactor**2 != D_factor:
            continue

        if (u_value + sqrt_Dfactor) % 2 != 0:
            continue

        p_candidate = (
            u_value + sqrt_Dfactor
        ) // 2

        q_candidate = (
            u_value - sqrt_Dfactor
        ) // 2

        if p_candidate * q_candidate == n:
            pair = tuple(
                sorted((p_candidate, q_candidate))
            )

            if pair not in result["factor_candidates"]:
                result["factor_candidates"].append(pair)

    if result["factor_candidates"]:
        result["success"] = True
        result["reason"] = "factorization recovered"
    else:
        result["reason"] = "no valid factor pair"

    return result


# ============================================================
# 6. BASIC TEST: 17 * 23
# ============================================================

print("\n" + "=" * 80)
print("5. BASIC FACTOR RECOVERY TEST")
print("=" * 80)

p0 = 17
q0 = 23
r0 = 5

n0 = p0 * q0
kappa30 = 1 - F(p0) * F(q0) * F(r0)

print("p =", p0)
print("q =", q0)
print("n =", n0)
print("r =", r0)
print("F(r) =", F(r0))
print("kappa3 =", kappa30)

basic_result = recover_factors(
    n0,
    kappa30,
    r0
)

print("\nRecovery:")
print(basic_result)

assert basic_result["success"]
assert (17, 23) in basic_result["factor_candidates"]

print("\nPASS: 17 * 23 recovered.")


# ============================================================
# 7. VERIFY THE p+q EQUATION DIRECTLY
# ============================================================

print("\n" + "=" * 80)
print("6. DIRECT p+q VERIFICATION")
print("=" * 80)

A0 = int(F(p0) * F(q0))
u0 = p0 + q0
v0 = p0 * q0

equation_value = (
    u0**2
    - (v0 + 1) * u0
    + (v0**2 - v0 + 1 - A0)
)

print("A =", A0)
print("u =", u0)
print("v =", v0)

print("\nQuadratic evaluated at actual p+q:")
print(equation_value)

assert equation_value == 0

print("\nPASS.")


# ============================================================
# 8. EXHAUSTIVE PRIME-PAIR TEST
# ============================================================

print("\n" + "=" * 80)
print("7. EXHAUSTIVE PRIME-PAIR RECOVERY")
print("=" * 80)

prime_limit = 100
primes = list(sp.primerange(2, prime_limit + 1))

total = 0
successes = 0
failures = []

for pp in primes:
    for qq in primes:

        if pp >= qq:
            continue

        n0 = pp * qq

        for rr in range(-10, 21):

            k3 = 1 - F(pp) * F(qq) * F(rr)

            result = recover_factors(
                n0,
                k3,
                rr
            )

            total += 1

            if result["success"]:
                pair = tuple(
                    sorted((pp, qq))
                )

                if pair in result["factor_candidates"]:
                    successes += 1
                else:
                    failures.append(
                        (pp, qq, rr, result)
                    )
            else:
                failures.append(
                    (pp, qq, rr, result)
                )

print("Prime limit:", prime_limit)
print("Total tests:", total)
print("Successful recoveries:", successes)
print("Failures:", len(failures))

if total:
    print(
        "Success rate:",
        successes / total
    )

if failures:
    print("\nFirst 10 failures:")
    for failure in failures[:10]:
        print(failure)
else:
    print("\nPASS: every tested instance recovered.")


# ============================================================
# 9. TEST CONTROLLED r VALUES
# ============================================================

print("\n" + "=" * 80)
print("8. CONTROLLED r VALUES")
print("=" * 80)

for rr in range(-5, 11):

    print(
        f"r={rr:3d}   "
        f"F(r)={F(rr):4d}"
    )


# ============================================================
# 10. SPECIAL VALUES r=0 AND r=1
# ============================================================

print("\n" + "=" * 80)
print("9. SPECIAL VALUES r=0 AND r=1")
print("=" * 80)

K3_r0 = sp.factor(K3.subs(r, 0))
K3_r1 = sp.factor(K3.subs(r, 1))

print("K3(r=0):")
print(K3_r0)

print("\nK3(r=1):")
print(K3_r1)

print("\nDifference:")
print(sp.factor(K3_r0 - K3_r1))

assert sp.expand(K3_r0 - K3_r1) == 0

print("\nPASS: r=0 and r=1 are identical.")


# ============================================================
# 11. SPECIAL CONTROL r=-1
# ============================================================

print("\n" + "=" * 80)
print("10. SPECIAL VALUE r=-1")
print("=" * 80)

print("F(-1) =", F(-1))

print("\nK3(r=-1):")
print(sp.factor(K3.subs(r, -1)))


# ============================================================
# 12. MODULAR STRUCTURE OF F(r)
# ============================================================

print("\n" + "=" * 80)
print("11. MODULAR ROOT STRUCTURE")
print("=" * 80)

for modulus in range(2, 31):

    roots = []

    for rr in range(modulus):

        if F(rr) % modulus == 0:
            roots.append(rr)

    print(
        f"mod {modulus:2d}: "
        f"roots of F(r)=0 -> {roots}"
    )


# ============================================================
# 13. MODULAR COLLAPSE SEARCH
# ============================================================

print("\n" + "=" * 80)
print("12. MODULAR COLLAPSE SEARCH")
print("=" * 80)

for modulus in range(2, 21):

    interesting = []

    for rr in range(modulus):

        values = set()

        for pp in range(modulus):
            for qq in range(modulus):

                k = (
                    1
                    - F(pp)
                    * F(qq)
                    * F(rr)
                ) % modulus

                values.add(k)

        # A very small image means strong collapse.
        if len(values) <= 3:

            interesting.append(
                (
                    rr,
                    F(rr) % modulus,
                    sorted(values)
                )
            )

    if interesting:

        print(
            "\nmod",
            modulus,
            "interesting r values:"
        )

        for item in interesting:
            print(
                "  r =",
                item[0],
                "F(r) =",
                item[1],
                "kappa residues =",
                item[2]
            )


# ============================================================
# 14. SEARCH FOR r WHERE F(r) HAS SMALL VALUES
# ============================================================

print("\n" + "=" * 80)
print("13. SMALL F(r) SEARCH")
print("=" * 80)

small_values = []

for rr in range(-100, 101):

    fr = F(rr)

    if abs(fr) <= 20:

        small_values.append(
            (rr, fr)
        )

print(
    "r values with |F(r)| <= 20:"
)

print(small_values)


# ============================================================
# 15. 3 -> 4 TRANSITION FOR SPECIAL s
# ============================================================

print("\n" + "=" * 80)
print("14. CONTROLLED s TRANSITION")
print("=" * 80)

for ss in [-2, -1, 0, 1, 2, 3]:

    transition = sp.factor(
        (K4 - K3).subs(s, ss)
    )

    print(
        f"\ns={ss}:"
    )

    print(transition)


# ============================================================
# 16. CHECK TRANSITION QUADRATIC IN s
# ============================================================

print("\n" + "=" * 80)
print("15. DEGREE IN CONTROLLED s")
print("=" * 80)

degree_s = sp.degree(
    sp.expand(K4 - K3),
    s
)

print("degree =", degree_s)

assert degree_s == 2

print("\nPASS: transition is quadratic in s.")


# ============================================================
# 17. COMPLETE SYMBOLIC FACTOR-RECOVERY FORMULA
# ============================================================

print("\n" + "=" * 80)
print("16. FINAL SYMBOLIC FACTOR-RECOVERY FORMULA")
print("=" * 80)

print(r"""
Define

    F(x) = x^2-x+1.

For

    n=pq

and controlled r,

    kappa3 = 1-F(p)F(q)F(r).

Therefore

    A = (1-kappa3)/F(r).

Let

    u=p+q.

Then

    A = u^2-(n+1)u+n^2-n+1.

Hence

    u^2-(n+1)u+n^2-n+1-A = 0.

The discriminant is

    D_u = 4A-3(n-1)^2.

Therefore

    p+q = [n+1 +/- sqrt(D_u)]/2.

Finally,

    p,q = [u +/- sqrt(u^2-4n)]/2.
""")


# ============================================================
# 18. FINAL AUTOMATED SUMMARY
# ============================================================

print("\n" + "=" * 80)
print("17. FINAL SUMMARY")
print("=" * 80)

print("""
EXACT IDENTITY
--------------

kappa3 = 1 - F(p)F(q)F(r)

where

F(x) = x^2-x+1.


CONTROLLED r
------------

A = (1-kappa3)/F(r)

gives

A = F(p)F(q).


SYMMETRIC REDUCTION
-------------------

A = u^2-(n+1)u+n^2-n+1

where

u=p+q
n=pq.


RECOVERY OF p+q
---------------

D_u = 4A - 3(n-1)^2

u = (n+1 +/- sqrt(D_u))/2.


RECOVERY OF FACTORS
-------------------

D_factor = u^2 - 4n

p,q = (u +/- sqrt(D_factor))/2.


CENTRAL RESEARCH QUESTION
-------------------------

Can kappa3, or an equivalent quantity A,
be obtained for a hidden factorization n=pq
without already knowing p and q?

If YES:
    this structure may provide a factor-recovery
    mechanism.

If NO:
    the result remains an exact algebraic
    characterization rather than a factoring
    algorithm.
""")

print("=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)