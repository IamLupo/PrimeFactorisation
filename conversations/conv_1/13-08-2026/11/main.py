#!/usr/bin/env python3
"""
==============================================================================
KAPPA NEXT EXPERIMENT
MULTI-AUXILIARY SECOND-INVARIANT / ELIMINATION SEARCH
==============================================================================

Purpose
-------
The previous experiments established:

    K_r = 1 - A*F(r)

Therefore all K_r values are evaluations of the same hidden scalar A.

This experiment deliberately searches for something beyond that trivial
one-dimensional structure.

We classify candidate relations as:

    TRIVIAL      : follows directly from K_r = 1 - A F(r)
    A-DEPENDENT  : still contains the hidden A
    N-ONLY       : computable from n alone
    U-DEPENDENT  : depends on the hidden u=p+q
    NEW          : survives elimination of A/K and produces a nontrivial
                   n/u relation

IMPORTANT
---------
The oracle values A and K are used ONLY for discovery/verification.
A candidate is NOT counted as an n-only result unless its final formula
contains no A, K, p, or q.

No CSV files are written.
All important results are printed.

Core identities
---------------
    F(x) = x^2 - x + 1
    n = p*q
    u = p+q

    A = F(p)F(q)

    A(u) = u^2 -(n+1)u + n^2-n+1

    K_r = 1 - A F(r)

We search:

  1. normalized K/F invariants
  2. cross-ratios of auxiliary values
  3. finite differences in F(r)
  4. products/determinants of multiple K values
  5. polynomial relations among (F(r), K_r)
  6. substitutions using F(r)=r²-r+1
  7. relations involving n and r
  8. elimination of A from symbolic-looking expressions
  9. multi-target recurrence/invariant tests
 10. search for a second independent scalar invariant

The most important output is the classification of identities.

==============================================================================
"""

from math import gcd, isqrt
from itertools import combinations, product
from fractions import Fraction
import random


# ============================================================================
# BASIC FUNCTIONS
# ============================================================================

def F(x):
    return x*x - x + 1


def bits(x):
    return abs(x).bit_length()


def factor_pair_from_primes(p, q):
    n = p*q
    u = p+q
    A = F(p)*F(q)
    return n, u, A


def A_from_n_u(n, u):
    return u*u - (n+1)*u + n*n - n + 1


def K_from_A(A, r):
    return 1 - A*F(r)


def discriminant_from_n_A(n, A):
    return 4*A - 3*(n-1)*(n-1)


# ============================================================================
# AUXILIARY SETS
# ============================================================================

AUX = {
    "first": (2, 3, 5, 7, 11, 13, 17, 19, 23),
    "odd_first": (3, 5, 7, 11, 13, 17, 19, 23),
    "every_other": (2, 5, 11, 17, 23, 31, 41, 47),
    "larger_first": (11, 13, 17, 19, 23, 29),
}


# ============================================================================
# TARGET GENERATION
# ============================================================================

def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d*d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def random_prime(bitsize):
    while True:
        x = random.getrandbits(bitsize)
        x |= (1 << (bitsize - 1))
        x |= 1
        if is_prime(x):
            return x


def generate_target(bitsize=25):
    p = random_prime(bitsize)
    q = random_prime(bitsize)

    while p == q:
        q = random_prime(bitsize)

    return factor_pair_from_primes(p, q) + (p, q)


# ============================================================================
# PRINT HEADER
# ============================================================================

print("=" * 78)
print("KAPPA NEXT EXPERIMENT")
print("MULTI-AUXILIARY SECOND-INVARIANT / ELIMINATION SEARCH")
print("=" * 78)

print("""
Core question:

    Do multiple auxiliary evaluations contain a second invariant
    beyond the single hidden scalar A?

We know:

    K_r = 1 - A F(r)

The experiment searches for relations that survive elimination of A
and are genuinely connected to n or u.
""")

print("AUXILIARY STRUCTURE")
print("-" * 78)

for name, rs in AUX.items():
    B = 1
    for r in rs:
        B *= F(r)

    print(
        f"{name:16s} count={len(rs):2d} "
        f"B_bits={bits(B):2d} B={B}"
    )
    print(f"  r={rs}")


# ============================================================================
# REPRESENTATIVE TARGET
# ============================================================================

random.seed(20260814)

n, u, A, p, q = generate_target(25)

print("\n" + "=" * 78)
print("REPRESENTATIVE TARGET")
print("=" * 78)

print(f"p       = {p}")
print(f"q       = {q}")
print(f"n       = {n}")
print(f"n bits  = {bits(n)}")
print(f"u       = {u}")
print(f"u bits  = {bits(u)}")
print(f"A       = {A}")
print(f"A bits  = {bits(A)}")


# ============================================================================
# SECTION 1: VERIFY BASE STRUCTURE
# ============================================================================

print("\n" + "=" * 78)
print("1. BASE AFFINE STRUCTURE")
print("=" * 78)

print("""
For every auxiliary r:

    K_r = 1 - A F(r)

Hence:

    (1-K_r)/F(r) = A

and

    (K_r-K_s)/(F(r)-F(s)) = -A.
""")

print(f"{'r':>4} {'F(r)':>14} {'K_r':>28} {'recovered A':>28}")

for r in AUX["first"]:
    fr = F(r)
    kr = K_from_A(A, r)
    recovered = Fraction(1-kr, fr)

    print(
        f"{r:4d} {fr:14d} {kr:28d} {str(recovered):>28}"
    )


# ============================================================================
# SECTION 2: SEARCH FOR SECOND SCALAR INVARIANTS
# ============================================================================

print("\n" + "=" * 78)
print("2. SEARCH FOR SECOND SCALAR INVARIANTS")
print("=" * 78)

rs = AUX["first"]

K = {r: K_from_A(A, r) for r in rs}
FR = {r: F(r) for r in rs}

# Candidate scalar quantities
candidate_scalars = {}

for r in rs:
    candidate_scalars[f"(1-K_{r})/F_{r}"] = Fraction(
        1-K[r], FR[r]
    )

for r, s in combinations(rs, 2):
    candidate_scalars[f"(K_{r}-K_{s})/(F_{r}-F_{s})"] = Fraction(
        K[r]-K[s],
        FR[r]-FR[s]
    )

# Check how many distinct values occur.
groups = {}

for name, value in candidate_scalars.items():
    groups.setdefault(value, []).append(name)

print(f"candidate scalar expressions = {len(candidate_scalars)}")
print(f"distinct scalar values         = {len(groups)}")

for value, names in sorted(groups.items(), key=lambda z: str(z[0])):
    print(f"\nvalue = {value}")
    for name in names[:8]:
        print(f"  {name}")

    if len(names) > 8:
        print(f"  ... {len(names)-8} more")


# ============================================================================
# SECTION 3: CROSS-RATIO TEST
# ============================================================================

print("\n" + "=" * 78)
print("3. CROSS-RATIO / PROJECTIVE INVARIANT TEST")
print("=" * 78)

print("""
If K is an affine function of F, then cross-ratios of K values
must equal cross-ratios of F values:

 CR(K1,K2;K3,K4) = CR(F1,F2;F3,F4)

This eliminates A completely.

The question is whether the resulting invariant contains n or u.
""")

cross_results = []

for r1, r2, r3, r4 in combinations(rs, 4):
    a = Fraction(K[r1]-K[r3], K[r2]-K[r3])
    b = Fraction(K[r1]-K[r4], K[r2]-K[r4])
    crK = a / b

    c = Fraction(FR[r1]-FR[r3], FR[r2]-FR[r3])
    d = Fraction(FR[r1]-FR[r4], FR[r2]-FR[r4])
    crF = c / d

    equal = (crK == crF)

    cross_results.append(equal)

    if len(cross_results) <= 12:
        print(
            f"({r1},{r2};{r3},{r4}) "
            f"CR_K={crK} "
            f"CR_F={crF} "
            f"equal={equal}"
        )

print(f"\nCross-ratio identities: "
      f"{sum(cross_results)}/{len(cross_results)}")


# ============================================================================
# SECTION 4: DETERMINANT SEARCH
# ============================================================================

print("\n" + "=" * 78)
print("4. MULTI-AUXILIARY DETERMINANT SEARCH")
print("=" * 78)

print("""
Because

    K_r = 1 - A F(r),

the points (F(r), K_r) lie on one affine line.

Therefore every 3x3 determinant

    | 1 F(r1) K(r1) |
    | 1 F(r2) K(r2) |
    | 1 F(r3) K(r3) |

must vanish.

We verify this and then test whether replacing F(r)
with expressions involving n or r creates anything new.
""")

det_zero = 0
det_total = 0

for r1, r2, r3 in combinations(rs, 3):
    mat = [
        [1, FR[r1], K[r1]],
        [1, FR[r2], K[r2]],
        [1, FR[r3], K[r3]],
    ]

    det = (
        mat[0][0] * (mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1])
        - mat[0][1] * (mat[1][0]*mat[2][2] - mat[1][2]*mat[2][0])
        + mat[0][2] * (mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0])
    )

    det_total += 1

    if det == 0:
        det_zero += 1

print(f"3x3 determinants tested = {det_total}")
print(f"zero determinants        = {det_zero}")


# ============================================================================
# SECTION 5: SEARCH EXPRESSIONS INVOLVING n
# ============================================================================

print("\n" + "=" * 78)
print("5. N-DEPENDENT ELIMINATION SEARCH")
print("=" * 78)

print("""
We now test whether multiplying/eliminating the K values with
simple n-dependent coefficients creates a quantity involving u.

Candidate coefficients include:

    n
    n+1
    n-1
    n^2-1
    F(r)
    r
    r-1
    r+1

The test is deliberately conservative:
a result is interesting only if it is not a constant determined
solely by the auxiliary r-values.
""")

n_coeffs = {
    "n": n,
    "n-1": n-1,
    "n+1": n+1,
    "n^2-1": n*n-1,
}

interesting = []

for r, s in combinations(rs, 2):

    for cname, c in n_coeffs.items():

        # Eliminate A from:
        #
        # K_r = 1-AF_r
        # K_s = 1-AF_s
        #
        # by forming:
        #
        # F_s K_r - F_r K_s = F_s-F_r.
        #
        # Then test n-scaled versions.

        lhs = c * (FR[s]*K[r] - FR[r]*K[s])
        rhs = c * (FR[s]-FR[r])

        residual = lhs-rhs

        if residual == 0:
            # This is algebraically trivial.
            classification = "TRIVIAL"
        else:
            classification = "NONZERO"

        if classification != "TRIVIAL":
            interesting.append(
                (r, s, cname, residual)
            )

print(f"n-dependent elimination tests = "
      f"{len(list(combinations(rs,2))) * len(n_coeffs)}")

print(f"nontrivial residuals = {len(interesting)}")

if interesting:
    for item in interesting[:20]:
        print(item)
else:
    print("No new n-dependent residual found.")


# ============================================================================
# SECTION 6: SEARCH FOR RELATIONS BETWEEN AUXILIARY INDEX r AND n
# ============================================================================

print("\n" + "=" * 78)
print("6. r / F(r) / n RELATION SEARCH")
print("=" * 78)

expressions = {}

for r in rs:
    fr = F(r)

    expressions.setdefault("r", []).append(r)
    expressions.setdefault("F(r)", []).append(fr)
    expressions.setdefault("F(r)-r", []).append(fr-r)
    expressions.setdefault("F(r)-(r^2)", []).append(fr-r*r)
    expressions.setdefault("F(r)+r", []).append(fr+r)
    expressions.setdefault("F(r) mod n", []).append(fr % n)

for name, vals in expressions.items():
    distinct = len(set(vals))

    print(
        f"{name:20s} distinct={distinct:2d} "
        f"values={vals[:9]}"
    )


# ============================================================================
# SECTION 7: u-POLYNOMIAL / AUXILIARY INTERACTION
# ============================================================================

print("\n" + "=" * 78)
print("7. u-POLYNOMIAL / AUXILIARY INTERACTION")
print("=" * 78)

print("""
We test whether A(u) has a special relationship with F(r)
that survives elimination of A.

For each r we examine:

    A(u) mod F(r)

and compare it against expressions formed only from n,u,r.

This is NOT considered an n-only result unless the hidden u
also disappears from the final relation.
""")

for r in rs:
    m = F(r)

    au_mod = A_from_n_u(n, u) % m

    candidates = {
        "u": u % m,
        "n": n % m,
        "u-r": (u-r) % m,
        "u+r": (u+r) % m,
        "n-u": (n-u) % m,
        "n+u": (n+u) % m,
        "u^2": (u*u) % m,
        "n*u": (n*u) % m,
        "n^2": (n*n) % m,
    }

    hits = [
        name for name, value in candidates.items()
        if value == au_mod
    ]

    print(
        f"r={r:2d} F={m:4d} "
        f"A(u) mod F={au_mod:4d} "
        f"matches={hits}"
    )


# ============================================================================
# SECTION 8: SECOND-INVARIANT SEARCH USING RATIOS
# ============================================================================

print("\n" + "=" * 78)
print("8. RATIO-OF-DIFFERENCES SEARCH")
print("=" * 78)

print("""
A first-order ratio gives only A:

    (K_r-K_s)/(F_r-F_s) = -A.

We therefore test ratios involving TWO independent differences:

    R = ((K_r-K_s)/(F_r-F_s))
        /
        ((K_t-K_v)/(F_t-F_v))

For the current model this must equal 1.

Any deviation would indicate a second degree of freedom.
""")

ratio_values = []

for r, s, t, v in combinations(rs, 4):
    d1 = Fraction(K[r]-K[s], FR[r]-FR[s])
    d2 = Fraction(K[t]-K[v], FR[t]-FR[v])

    ratio = d1/d2
    ratio_values.append(ratio)

print(f"ratios tested = {len(ratio_values)}")
print(f"distinct ratios = {len(set(ratio_values))}")

print("sample ratios:")
for value in ratio_values[:15]:
    print(f"  {value}")


# ============================================================================
# SECTION 9: NONLINEAR K COMBINATIONS
# ============================================================================

print("\n" + "=" * 78)
print("9. NONLINEAR K-COMBINATION SEARCH")
print("=" * 78)

print("""
Now test low-degree combinations:

    K_r*K_s
    K_r + K_s
    K_r*K_s - K_t
    (1-K_r)(1-K_s)
    (K_r-1)(K_s-1)

and ask whether A can be eliminated between several such
expressions.

This is where a genuinely nonlinear second invariant would
first have a chance to appear.
""")

nonlinear_hits = []

for r, s in combinations(rs, 2):
    for t in rs:

        if t in (r, s):
            continue

        x = (K[r]-1)*(K[s]-1)
        y = (K[t]-1)

        # Since K_r-1 = -A F(r),
        # x / y = -A F(r)F(s)/F(t).
        #
        # Normalize the auxiliary factor away.

        normalized = Fraction(
            x * FR[t],
            y * FR[r] * FR[s]
        )

        if normalized == -A:
            nonlinear_hits.append(
                (r, s, t, normalized)
            )

print(f"nonlinear normalized identities = {len(nonlinear_hits)}")

for item in nonlinear_hits[:15]:
    print(
        f"r,s,t=({item[0]},{item[1]},{item[2]}) "
        f"value={item[3]}"
    )


# ============================================================================
# SECTION 10: MULTI-TARGET INDEPENDENCE TEST
# ============================================================================

print("\n" + "=" * 78)
print("10. MULTI-TARGET SECOND-INVARIANT TEST")
print("=" * 78)

print("""
A real structural invariant should survive changing p and q.

We generate several independent targets and compare normalized
auxiliary expressions.

If a quantity depends only on the auxiliary sequence, it should
be identical across targets.

If it depends on A, it should vary with A.

If it depends on n/u in a new way, its behavior should correlate
with n/u rather than simply A.
""")

targets = []

for i in range(8):
    target = generate_target(18)
    tn, tu, tA, tp, tq = target

    targets.append({
        "p": tp,
        "q": tq,
        "n": tn,
        "u": tu,
        "A": tA,
    })

    print(
        f"target {i+1}: "
        f"n_bits={bits(tn):2d} "
        f"u_bits={bits(tu):2d} "
        f"A_bits={bits(tA):2d}"
    )


# Candidate normalized quantities
def target_signature(target):
    ta = target["A"]

    out = {}

    for r, s in combinations(rs[:5], 2):
        kr = K_from_A(ta, r)
        ks = K_from_A(ta, s)

        value = Fraction(
            kr-ks,
            F(r)-F(s)
        )

        out[f"slope_{r}_{s}"] = value

    return out


signatures = [target_signature(t) for t in targets]

print("\nSlope behavior:")
for key in list(signatures[0])[:5]:
    vals = [s[key] for s in signatures]

    print(
        f"{key:20s} "
        f"distinct={len(set(vals))} "
        f"values={vals[:4]}"
    )


# ============================================================================
# SECTION 11: SEARCH FOR A SECOND DEGREE OF FREEDOM
# ============================================================================

print("\n" + "=" * 78)
print("11. DEGREE-OF-FREEDOM DIAGNOSTIC")
print("=" * 78)

print("""
For each target we treat the vector

    (K_r)

as a point in auxiliary-value space.

If all vectors lie on a one-dimensional affine family

    K = 1 - A F,

then changing p,q only changes one scalar coordinate A.

We numerically verify this by checking whether every pairwise
difference ratio is exactly determined by F(r).
""")

dimension_failures = 0
dimension_tests = 0

for target in targets:
    ta = target["A"]

    kt = {r: K_from_A(ta, r) for r in rs}

    for r, s, t in combinations(rs[:6], 3):

        lhs = Fraction(
            kt[r]-kt[s],
            kt[r]-kt[t]
        )

        rhs = Fraction(
            F(r)-F(s),
            F(r)-F(t)
        )

        dimension_tests += 1

        if lhs != rhs:
            dimension_failures += 1

print(f"dimension tests = {dimension_tests}")
print(f"failures        = {dimension_failures}")

if dimension_failures == 0:
    print("PASS: all auxiliary vectors remain one-dimensional.")


# ============================================================================
# SECTION 12: SEARCH FOR n/u CORRELATION AFTER REMOVING A
# ============================================================================

print("\n" + "=" * 78)
print("12. SEARCH FOR n/u INFORMATION AFTER REMOVING A")
print("=" * 78)

print("""
We construct quantities where the entire A-dependent component
has been removed.

Examples:

    slope differences
    cross-ratio differences
    normalized products
    determinant residuals

Then we test whether any surviving quantity varies with n or u.

If every surviving quantity is constant across targets, the
auxiliary family has no detectable second invariant.
""")

survivors = []

for target in targets:

    tn = target["n"]
    tu = target["u"]
    ta = target["A"]

    kt = {r: K_from_A(ta, r) for r in rs}

    # Cross-ratio
    r1, r2, r3, r4 = rs[:4]

    cr = (
        Fraction(kt[r1]-kt[r3], kt[r2]-kt[r3])
        /
        Fraction(kt[r1]-kt[r4], kt[r2]-kt[r4])
    )

    cr_expected = (
        Fraction(F(r1)-F(r3), F(r2)-F(r3))
        /
        Fraction(F(r1)-F(r4), F(r2)-F(r4))
    )

    survivors.append({
        "n": tn,
        "u": tu,
        "cross_residual": cr-cr_expected,
    })

for i, row in enumerate(survivors, 1):
    print(
        f"target {i}: "
        f"cross-ratio residual={row['cross_residual']}"
    )


# ============================================================================
# SECTION 13: SEARCH FOR POLYNOMIAL RELATIONS IN n,u,r
# ============================================================================

print("\n" + "=" * 78)
print("13. SMALL POLYNOMIAL RELATION SEARCH")
print("=" * 78)

print("""
Search:

    c0 + c1*r + c2*u + c3*n
      + c4*r*u + c5*r*n
      + c6*u^2 + c7*n*u + c8*n^2

for simple repeated zero relations.

A relation is interesting only if:

    * it vanishes across many targets,
    * it involves n/u/r,
    * and it is not an obvious defining identity.
""")

# Candidate expressions.
def polynomial_candidates(nv, uv, r):
    return {
        "1": 1,
        "r": r,
        "u": uv,
        "n": nv,
        "r*u": r*uv,
        "r*n": r*nv,
        "u^2": uv*uv,
        "n*u": nv*uv,
        "n^2": nv*nv,
        "r^2": r*r,
        "r*u^2": r*uv*uv,
        "r*n*u": r*nv*uv,
    }


zero_counts = {}

for target in targets:
    tn = target["n"]
    tu = target["u"]

    for r in rs:
        vals = polynomial_candidates(tn, tu, r)

        # Test pairwise differences against zero.
        # This intentionally looks only for direct zero identities.
        for name, value in vals.items():
            zero_counts[name] = zero_counts.get(name, 0) + (value == 0)

print("Direct zero counts:")
for name, count in sorted(zero_counts.items(), key=lambda x: -x[1]):
    if count:
        print(
            f"{name:12s} {count}/{len(targets)*len(rs)}"
        )

print("\nNo nontrivial direct-zero invariant is expected here.")


# ============================================================================
# SECTION 14: FINAL CLASSIFICATION
# ============================================================================

print("\n" + "=" * 78)
print("FINAL CLASSIFICATION")
print("=" * 78)

print("""
RESULT CATEGORIES
-----------------

A. AFFINE / TRIVIAL
-------------------
    K_r = 1-AF(r)

    All pairwise slopes equal -A.
    All 3x3 determinants vanish.
    All cross-ratios match those of F(r).

B. A-ELIMINATION
-----------------
    F(s)K_r - F(r)K_s = F(s)-F(r)

This removes A but leaves only an auxiliary identity.

C. SECOND-INVARIANT SEARCH
---------------------------
The experiment searched for a quantity that:

    * survives elimination of A,
    * is not determined solely by r,
    * depends on n/u,
    * persists across targets.

D. N-ONLY ROUTE
----------------
A result counts as n-only ONLY if the final expression can be
computed without p, q, u, A, or K.

E. IMPORTANT NEGATIVE RESULT
----------------------------
If the auxiliary vector remains exactly one-dimensional across
all targets, then adding more auxiliary primes cannot create
a second independent scalar by itself.

The only possible breakthrough would therefore have to come from
an additional relation connecting the auxiliary construction to n,
rather than from K-elimination alone.
""")

# ============================================================================
# SUMMARY NUMBERS
# ============================================================================

print("=" * 78)
print("SUMMARY")
print("=" * 78)

print(f"Representative n bits : {bits(n)}")
print(f"Representative u bits : {bits(u)}")
print(f"Representative A bits : {bits(A)}")

print(
    f"Cross-ratio identities : "
    f"{sum(cross_results)}/{len(cross_results)}"
)

print(
    f"3x3 determinant zeros  : "
    f"{det_zero}/{det_total}"
)

print(
    f"Second-dimension tests : "
    f"{dimension_tests}"
)

print(
    f"Second-dimension failures : "
    f"{dimension_failures}"
)

if dimension_failures == 0:
    print("""
FINAL RESULT:

The tested auxiliary family remains exactly one-dimensional:

    K_vector = 1 - A * F_vector

No second independent auxiliary invariant was detected.

The next useful direction is NOT more auxiliary primes.
It is searching for an independent relation connecting this
one-dimensional auxiliary family to n.
""")
else:
    print("""
POTENTIALLY INTERESTING:

A deviation from the one-dimensional affine model was detected.
Inspect the reported residuals before proceeding.
""")

print("=" * 78)
print("DONE")
print("=" * 78)
