#!/usr/bin/env python3

"""
==============================================================================
KAPPA NEXT EXPERIMENT
DISCRIMINANT / SUM-OF-FACTORS / MULTI-MODULUS S-RECOVERY SEARCH
==============================================================================

Core question:

    Can n alone, together with the auxiliary moduli

        F(r) = r^2-r+1

    constrain

        s = p+q

    through the discriminant

        Delta = s^2 - 4n = (p-q)^2 ?

We deliberately DO NOT use p or q in the discovery stage.

The experiment has four layers:

    1. Study s^2 modulo F(r).
    2. Determine whether n modulo F(r) restricts possible s residues.
    3. Combine the restrictions over many auxiliary moduli.
    4. Test whether the resulting s candidates correlate with the true
       factor sum across independent targets.

No CSV.
Only printed output.
==============================================================================

"""

from math import gcd, isqrt
from itertools import combinations, product
from collections import defaultdict

SEED = 20260814

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17,
    19, 23, 29, 31, 37, 41, 43, 47
]

# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------

TARGETS = [
    (735067, 956861),
    (21390503, 22630357),
    (17543153, 31300123),
    (22612043, 26706517),

    (563, 587),
    (719, 761),
    (1009, 1031),
    (1423, 1451),
    (2003, 2017),
    (3001, 3011),
    (4001, 4013),
    (5003, 5021),
]

# ---------------------------------------------------------------------------
# Basic functions
# ---------------------------------------------------------------------------

def F(r):
    return r * r - r + 1


def sixth_root_residues(m):
    return [x for x in range(m) if pow(x, 6, m) == 1]


def quadratic_roots_mod(n, m):
    """
    Return x mod m satisfying x^2 == n mod m.

    Brute force is intentional here because auxiliary moduli are small.
    """
    return [x for x in range(m) if (x * x - n) % m == 0]


def canonical_pair(a, m):
    """
    Normalize x and -x as one discriminant-sign pair.
    """
    b = (-a) % m
    return min(a, b), max(a, b)


def residue_signature(value, moduli):
    return tuple(value % m for m in moduli)


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------

print("=" * 78)
print("KAPPA NEXT EXPERIMENT")
print("DISCRIMINANT / SUM-OF-FACTORS / MULTI-MODULUS S-RECOVERY SEARCH")
print("=" * 78)

print()
print(f"random seed = {SEED}")
print(f"R values    = {R_VALUES}")
print(f"targets     = {len(TARGETS)}")

# ---------------------------------------------------------------------------
# 1. Representative target
# ---------------------------------------------------------------------------

p0, q0 = TARGETS[0]
n0 = p0 * q0
s0 = p0 + q0
delta0 = (p0 - q0) ** 2

print()
print("=" * 78)
print("1. REPRESENTATIVE TARGET")
print("=" * 78)

print(f"p       = {p0}")
print(f"q       = {q0}")
print(f"n       = {n0}")
print(f"s=p+q   = {s0}")
print(f"Delta   = {delta0}")
print(f"n bits  = {n0.bit_length()}")
print(f"s bits  = {s0.bit_length()}")

# ---------------------------------------------------------------------------
# 2. Auxiliary moduli
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("2. AUXILIARY MODULI")
print("=" * 78)

MODULI = []

for r in R_VALUES:
    m = F(r)
    MODULI.append(m)

    roots = sixth_root_residues(m)

    print(
        f"r={r:2d} "
        f"F={m:6d} "
        f"factor={m:4d} "
        f"sixth_roots={len(roots):3d}"
    )

# ---------------------------------------------------------------------------
# 3. Fundamental discriminant relation
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("3. DISCRIMINANT RELATION")
print("=" * 78)

print(
    """
We have

    s = p+q
    n = pq

and therefore

    s^2 - 4n = (p-q)^2.

Modulo m=F(r):

    s^2 = 4n + Delta (mod m).

The key question is whether n alone gives enough information to
restrict s or Delta modulo the auxiliary moduli.
"""
)

for r in R_VALUES[:9]:
    m = F(r)

    nmod = n0 % m
    smod = s0 % m
    dmod = delta0 % m

    print(
        f"r={r:2d} m={m:6d} "
        f"n={nmod:6d} "
        f"s={smod:6d} "
        f"Delta={dmod:6d} "
        f"check={(smod*smod - 4*nmod) % m:6d}"
    )

# ---------------------------------------------------------------------------
# 4. N-ONLY possible s residues
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("4. N-ONLY S-RESIDUE SEARCH")
print("=" * 78)

print(
    """
For each modulus m we enumerate all s residues.

Without additional information, every s mod m is possible.

We then impose only the necessary condition that

    Delta = s^2 - 4n

must be a quadratic residue modulo m.

This is a deliberately weak condition, because Delta must globally be
the square (p-q)^2.

If this filtering becomes strong when several moduli are combined,
that is worth pursuing.
"""
)

for r in R_VALUES[:9]:
    m = F(r)
    nmod = n0 % m

    valid_s = []
    valid_delta = []

    square_set = {
        (x * x) % m
        for x in range(m)
    }

    for s in range(m):
        d = (s * s - 4 * nmod) % m
        if d in square_set:
            valid_s.append(s)
            valid_delta.append(d)

    print(
        f"r={r:2d} m={m:6d} "
        f"n_mod={nmod:6d} "
        f"valid_s={len(valid_s):6d}/{m:6d} "
        f"valid_Delta={len(set(valid_delta)):6d}"
    )

# ---------------------------------------------------------------------------
# 5. Stronger condition: s and Delta must correspond to a factor pair
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("5. MODULAR FACTOR-PAIR CONSISTENCY")
print("=" * 78)

print(
    """
Modulo m, a genuine factor pair satisfies

    x*y = n
    x+y  = s.

Therefore x and y are roots of

    X^2 - sX + n = 0.

The discriminant condition is exactly

    s^2 - 4n is a square mod m.

We enumerate the possible s residues and record the number of
corresponding root pairs.
"""
)

for r in R_VALUES[:9]:
    m = F(r)
    nmod = n0 % m

    valid = []

    for s in range(m):
        roots = [
            x for x in range(m)
            if (x * x - s * x + nmod) % m == 0
        ]

        if roots:
            valid.append((s, roots))

    print(
        f"r={r:2d} m={m:6d} "
        f"possible_s={len(valid):6d}"
    )

    if m <= 50:
        print(f"  s candidates = {[x[0] for x in valid]}")

# ---------------------------------------------------------------------------
# 6. Compare actual s against n-only candidates
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("6. TRUE-S MEMBERSHIP CONTROL")
print("=" * 78)

for r in R_VALUES[:9]:
    m = F(r)
    nmod = n0 % m
    true_s = s0 % m

    possible = []

    for s in range(m):
        if any(
            (x * x - s * x + nmod) % m == 0
            for x in range(m)
        ):
            possible.append(s)

    print(
        f"r={r:2d} m={m:6d} "
        f"true_s={true_s:6d} "
        f"in_candidates={true_s in possible} "
        f"candidate_count={len(possible):6d}"
    )

# ---------------------------------------------------------------------------
# 7. Multi-modulus signature search
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("7. MULTI-MODULUS S-SIGNATURE COLLAPSE")
print("=" * 78)

print(
    """
Now combine several auxiliary moduli.

For a chosen subset M=(m1,...,mk), we compute all residue vectors

    (s mod m1, ..., s mod mk)

that satisfy the modular factor-pair condition for every modulus.

The goal is not to reconstruct the exact integer s yet.

The goal is to see whether the admissible residue signatures become
substantially smaller than the naive product of the moduli.
"""
)

for k in [1, 2, 3, 4, 5]:
    chosen = MODULI[:k]

    counts = []

    for target_index, (p, q) in enumerate(TARGETS):
        n = p * q

        candidate_sets = []

        for m in chosen:
            nmod = n % m

            vals = []

            for s in range(m):
                ok = False

                for x in range(m):
                    y = (s - x) % m

                    if (x * y - nmod) % m == 0:
                        ok = True
                        break

                if ok:
                    vals.append(s)

            candidate_sets.append(vals)

        count = 1
        for vals in candidate_sets:
            count *= len(vals)

        counts.append(count)

    print(
        f"k={k} "
        f"moduli={chosen} "
        f"candidate_signature_counts="
        f"min={min(counts)} "
        f"max={max(counts)} "
        f"avg={sum(counts)/len(counts):.1f}"
    )

# ---------------------------------------------------------------------------
# 8. CRT reconstruction of candidate s residues
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("8. CRT CANDIDATE-S INTERVAL TEST")
print("=" * 78)

print(
    """
For pairwise-coprime auxiliary moduli, CRT would combine the s residues.

The F(r) are not all coprime, so we use generalized CRT consistency.

Instead of enumerating all huge integers, we construct compatible
residue classes incrementally.
"""
)


def crt_pair(a1, m1, a2, m2):
    """
    Generalized CRT.

    Returns (a,m) with
        x == a1 mod m1
        x == a2 mod m2

    or None if inconsistent.
    """

    g = gcd(m1, m2)

    if (a2 - a1) % g != 0:
        return None

    m1g = m1 // g
    m2g = m2 // g

    if m2g == 1:
        t = 0
    else:
        inv = pow(m1g, -1, m2g)
        t = ((a2 - a1) // g * inv) % m2g

    a = a1 + m1 * t
    m = m1 * m2g

    return a % m, m


def candidate_s_residues(n, moduli):
    """
    Build generalized-CRT compatible s classes satisfying

        x*y == n mod m

    for each auxiliary modulus.
    """

    classes = [(0, 1)]

    for m in moduli:

        nmod = n % m

        possible_s = []

        for s in range(m):

            found = False

            for x in range(m):
                y = (s - x) % m

                if (x * y - nmod) % m == 0:
                    found = True
                    break

            if found:
                possible_s.append(s)

        new_classes = []

        for a, mod in classes:
            for s in possible_s:

                merged = crt_pair(a, mod, s, m)

                if merged is not None:
                    new_classes.append(merged)

        # Deduplicate.
        classes = list(set(new_classes))

        # Keep experiment bounded.
        if len(classes) > 200000:
            classes = classes[:200000]

    return classes


for k in [2, 3, 4, 5]:
    chosen = MODULI[:k]

    classes = candidate_s_residues(n0, chosen)

    modulus = 1

    for m in chosen:
        modulus = modulus * m // gcd(modulus, m)

    print(
        f"k={k} "
        f"moduli={chosen} "
        f"CRT_classes={len(classes)} "
        f"combined_modulus={modulus} "
        f"modulus_bits={modulus.bit_length()}"
    )

# ---------------------------------------------------------------------------
# 9. Does n determine s modulo the combined auxiliary modulus?
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("9. N-ONLY COLLISION TEST FOR S")
print("=" * 78)

print(
    """
Across targets we group by the complete n-residue vector

    (n mod F(r1), ..., n mod F(rk))

and ask how many different true s-residue vectors occur.

If the ambiguity collapses to one, then n determines s modulo the
combined auxiliary modulus for that k.
"""
)

for k in [1, 2, 3, 4, 5]:

    chosen = MODULI[:k]

    groups = defaultdict(set)

    for p, q in TARGETS:

        n = p * q
        s = p + q

        n_sig = residue_signature(n, chosen)
        s_sig = residue_signature(s, chosen)

        groups[n_sig].add(s_sig)

    ambiguous = [
        (key, values)
        for key, values in groups.items()
        if len(values) > 1
    ]

    print(
        f"k={k} "
        f"groups={len(groups):3d} "
        f"ambiguous={len(ambiguous):3d}"
    )

    if ambiguous:
        print(
            "  sample="
            + str(ambiguous[:3])
        )

# ---------------------------------------------------------------------------
# 10. Discriminant collision test
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("10. N-ONLY COLLISION TEST FOR DISCRIMINANT")
print("=" * 78)

print(
    """
The same test is performed for

    Delta = (p-q)^2.

If equal n-signatures repeatedly imply equal Delta-signatures,
that would be much stronger than simply studying s.
"""
)

for k in [1, 2, 3, 4, 5]:

    chosen = MODULI[:k]

    groups = defaultdict(set)

    for p, q in TARGETS:

        n = p * q
        delta = (p - q) ** 2

        n_sig = residue_signature(n, chosen)
        d_sig = residue_signature(delta, chosen)

        groups[n_sig].add(d_sig)

    ambiguous = [
        (key, values)
        for key, values in groups.items()
        if len(values) > 1
    ]

    print(
        f"k={k} "
        f"groups={len(groups):3d} "
        f"ambiguous_Delta={len(ambiguous):3d}"
    )

# ---------------------------------------------------------------------------
# 11. Search for direct n -> s polynomial maps
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("11. LOW-DEGREE N -> S MODULAR MAP SEARCH")
print("=" * 78)

print(
    """
Search whether s modulo F(r) can be represented by a low-degree
polynomial in n modulo F(r):

    s = c0 + c1*n + c2*n^2 + ...

This is only a diagnostic. A successful polynomial that persists
across targets would be a major signal.
"""
)

for degree in range(1, 5):

    print()
    print(f"degree={degree}")

    for r in R_VALUES[:9]:

        m = F(r)

        data = []

        for p, q in TARGETS:
            n = p * q
            s = p + q

            data.append(
                (n % m, s % m)
            )

        # Group equal n residues.
        grouped = defaultdict(set)

        for x, y in data:
            grouped[x].add(y)

        # A deterministic function s=f(n) exists iff every n residue
        # observed has exactly one s residue.
        deterministic = all(
            len(vals) == 1
            for vals in grouped.values()
        )

        print(
            f"r={r:2d} m={m:6d} "
            f"deterministic={deterministic} "
            f"n_classes={len(grouped):3d}"
        )

# ---------------------------------------------------------------------------
# 12. Search polynomial discriminant identities
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("12. DISCRIMINANT POLYNOMIAL IDENTITY SEARCH")
print("=" * 78)

print(
    """
Test simple candidates for a relation

    Delta = P(n)

modulo F(r).

Candidates include:
    n
    n-1
    n+1
    n^2
    n^2-1
    n^2-n
    n^2+n
    n^2-n+1
    n^2+n+1
"""
)

def candidate_polynomials(n):
    return {
        "n": n,
        "n-1": n - 1,
        "n+1": n + 1,
        "n^2": n * n,
        "n^2-1": n * n - 1,
        "n^2-n": n * n - n,
        "n^2+n": n * n + n,
        "n^2-n+1": n * n - n + 1,
        "n^2+n+1": n * n + n + 1,
    }


for name in candidate_polynomials(n0):

    hits = 0
    total = 0

    for p, q in TARGETS:

        n = p * q
        delta = (p - q) ** 2

        for m in MODULI[:9]:

            total += 1

            if delta % m == candidate_polynomials(n)[name] % m:
                hits += 1

    print(
        f"{name:12s} "
        f"hits={hits:4d}/{total:4d} "
        f"fraction={hits/total:.4f}"
    )

# ---------------------------------------------------------------------------
# 13. Factor-side F(p)+F(q) as an s detector
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("13. CYCLOTOMIC SUM AS AN s-DETECTOR")
print("=" * 78)

print(
    """
For F(x)=x^2-x+1:

    F(p)+F(q)
      = s^2 - s - 2n + 2.

Therefore

    s^2-s
      = F(p)+F(q)+2n-2.

This experiment asks whether auxiliary information can determine
the left side from n alone.

We compare the factor-side quantity against n-only candidates.
"""
)

for name, func in [
    ("n", lambda n: n),
    ("n^2", lambda n: n*n),
    ("n^2-n", lambda n: n*n-n),
    ("n^2+n", lambda n: n*n+n),
    ("n^2-n+1", lambda n: n*n-n+1),
    ("n^2+n+1", lambda n: n*n+n+1),
]:

    hits = 0
    total = 0

    for p, q in TARGETS:

        n = p * q
        s = p + q

        lhs = s * s - s
        rhs = func(n)

        for m in MODULI[:9]:

            total += 1

            if lhs % m == rhs % m:
                hits += 1

    print(
        f"{name:12s} "
        f"hits={hits:4d}/{total:4d} "
        f"fraction={hits/total:.4f}"
    )

# ---------------------------------------------------------------------------
# 14. Resultant-style test
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("14. RESULTANT-STYLE ELIMINATION")
print("=" * 78)

print(
    """
The factor polynomial is

    X^2 - sX + n.

The auxiliary cyclotomic polynomial is

    X^2 - X + 1.

A common factor root would imply a resultant relation between n and s.

We compute the resultant symbolically by direct algebra.

For

    f(X)=X^2-sX+n
    g(X)=X^2-X+1

the resultant is:

    Res_X(f,g)
      = n^2 + n(1-2s) + (s^2-s+1).

The question:

    Does this resultant vanish or factor modulo F(r)
    in a target-independent way?
"""
)

def resultant_value(n, s):
    return (
        n*n
        + n*(1 - 2*s)
        + (s*s - s + 1)
    )


zero_counts = []

for r in R_VALUES[:9]:

    m = F(r)

    hits = 0

    for p, q in TARGETS:

        n = p*q
        s = p+q

        if resultant_value(n, s) % m == 0:
            hits += 1

    zero_counts.append(hits)

    print(
        f"r={r:2d} m={m:6d} "
        f"resultant_zero_targets={hits}/{len(TARGETS)}"
    )

# ---------------------------------------------------------------------------
# 15. Search whether resultant zero predicts factor alignment
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("15. RESULTANT ROOT-ALIGNMENT CONTROL")
print("=" * 78)

print(
    """
When the resultant is zero modulo m, the two quadratics share a root.

That means some residue x satisfies simultaneously:

    x^2-sx+n = 0
    x^2-x+1 = 0.

Equivalently, a factor residue may lie on the cyclotomic curve.

We test whether resultant-zero events occur substantially more often
than expected by chance.
"""
)

for r in R_VALUES[:9]:

    m = F(r)

    result_hits = 0
    direct_hits = 0

    for p, q in TARGETS:

        n = p*q
        s = p+q

        if resultant_value(n, s) % m == 0:
            result_hits += 1

        roots = [
            x for x in range(m)
            if (x*x - x + 1) % m == 0
            and (x*x - s*x + n) % m == 0
        ]

        if roots:
            direct_hits += 1

    print(
        f"r={r:2d} m={m:6d} "
        f"resultant={result_hits:2d} "
        f"direct_common_root={direct_hits:2d}"
    )

# ---------------------------------------------------------------------------
# 16. Multi-target s-signature entropy
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("16. MULTI-TARGET S-SIGNATURE INFORMATION")
print("=" * 78)

print(
    """
For each k, compare:

    H(n-signature)

against

    H(s-signature | n-signature).

We use simple exact ambiguity counts rather than statistical entropy.

A value of 1 means the observed n-signature uniquely determines
the observed s-signature across the target set.
"""
)

for k in range(1, 10):

    chosen = MODULI[:k]

    groups = defaultdict(set)

    for p, q in TARGETS:

        n = p*q
        s = p+q

        nsig = residue_signature(n, chosen)
        ssig = residue_signature(s, chosen)

        groups[nsig].add(ssig)

    total_ambiguity = sum(
        len(v)
        for v in groups.values()
    )

    max_ambiguity = max(
        len(v)
        for v in groups.values()
    )

    print(
        f"k={k:2d} "
        f"groups={len(groups):3d} "
        f"total_s_signatures={total_ambiguity:3d} "
        f"max_ambiguity={max_ambiguity:2d}"
    )

# ---------------------------------------------------------------------------
# 17. Direct comparison with true factor sum
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("17. REPRESENTATIVE S RECOVERY CONTROL")
print("=" * 78)

for k in [1, 2, 3, 4, 5]:

    chosen = MODULI[:k]

    classes = candidate_s_residues(n0, chosen)

    combined_modulus = 1

    for m in chosen:
        combined_modulus = (
            combined_modulus * m
            // gcd(combined_modulus, m)
        )

    true_residue = s0 % combined_modulus

    contains = any(
        a == true_residue
        for a, mod in classes
        if mod == combined_modulus
    )

    print(
        f"k={k} "
        f"combined_modulus={combined_modulus} "
        f"true_s_residue={true_residue} "
        f"classes={len(classes)} "
        f"true_s_present={contains}"
    )

# ---------------------------------------------------------------------------
# 18. Final classification
# ---------------------------------------------------------------------------

print()
print("=" * 78)
print("18. FINAL CLASSIFICATION")
print("=" * 78)

print(
    """
A. DISCRIMINANT ROUTE
---------------------
    Delta = s^2 - 4n

The experiment tests whether auxiliary moduli restrict Delta or s
using n alone.

B. SUM ROUTE
------------
    s = p+q

If n-residue signatures determine s-residue signatures, then the
missing symmetric variable may become accessible modulo a large
combined modulus.

C. RESULTANT ROUTE
------------------
    Res_X(X^2-sX+n, X^2-X+1)

A repeated zero modulo F(r) would indicate that the factor polynomial
and auxiliary cyclotomic polynomial share a residue root.

D. BREAKTHROUGH CONDITION
-------------------------
A result becomes genuinely interesting only if:

    1. It uses n and auxiliary F(r).
    2. It does not use p, q, s, Delta, u, A, or K to discover it.
    3. It survives across independent targets.
    4. It reduces the possible s/Delta residue classes substantially.
    5. The reduction is stronger than a generic quadratic-residue
       or CRT counting effect.

E. IMPORTANT CONTROL
--------------------
Even if s is determined modulo some auxiliary modulus, that does NOT
yet factor n.

The next threshold would be:

    combined_modulus > plausible range for s

or a sufficiently strong congruence that leaves only a tiny number
of integer s candidates.

Then:

    Delta = s^2 - 4n

can be tested for being a perfect square.

That would turn the modular relation into an actual factorization
mechanism.
"""
)

print()
print("Targets tested :", len(TARGETS))
print("Moduli tested  :", len(MODULI))
print("CSV output     : NONE")
print("=" * 78)
print("DONE")
print("=" * 78)

