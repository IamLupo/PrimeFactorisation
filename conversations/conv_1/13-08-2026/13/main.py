"""
==============================================================================
KAPPA NEXT EXPERIMENT
SYMMETRIC FACTOR / CYCLOTOMIC PRODUCT / n-ONLY ELIMINATION SEARCH
==============================================================================

Goal:

Previous experiments showed that

    F(r) = r^2-r+1

has the intrinsic sixth-root property

    r^6 = 1 (mod F(r))

but n alone did not determine A or u modulo F(r).

This experiment changes direction.

Instead of searching arbitrary polynomials in n, we use the fact that

    n = p*q

and investigate whether the symmetric factor structure

    p+q, pq=n

creates a hidden relation involving

    F(p), F(q), F(r)

which can be eliminated down to an expression involving n and r alone.

The experiment has two modes:

    1. DIAGNOSTIC:
       p,q are known and we inspect factor-side quantities.

    2. ELIMINATION:
       eliminate s=p+q algebraically and search for relations
       depending only on n and r.

No CSV files.
Everything is printed.
==============================================================================

"""

from math import gcd
from itertools import combinations, product
from fractions import Fraction
import random


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814
random.seed(SEED)

R_VALUES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

TARGET_COUNTS = {
    20: 10,
    30: 10,
    40: 10,
    50: 10,
}

PRINT_LIMIT = 20


# ============================================================================
# BASIC ARITHMETIC
# ============================================================================

def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def random_prime(bits):
    while True:
        x = random.getrandbits(bits)
        x |= (1 << (bits - 1))
        x |= 1
        if is_prime(x):
            return x


def F(x):
    return x * x - x + 1


def poly_values(x):
    return {
        "x": x,
        "x-1": x - 1,
        "x+1": x + 1,
        "x^2-1": x * x - 1,
        "F(x)": F(x),
        "x^2+x+1": x * x + x + 1,
    }


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets():
    targets = []

    for bits, count in TARGET_COUNTS.items():
        for _ in range(count):
            pb = bits // 2
            qb = bits - pb

            p = random_prime(pb)
            q = random_prime(qb)

            if p == q:
                continue

            if p > q:
                p, q = q, p

            n = p * q
            s = p + q

            targets.append({
                "p": p,
                "q": q,
                "n": n,
                "s": s,
                "bits": n.bit_length(),
            })

    return targets


targets = generate_targets()


# ============================================================================
# HEADER
# ============================================================================

print("=" * 78)
print("KAPPA NEXT EXPERIMENT")
print("SYMMETRIC FACTOR / CYCLOTOMIC PRODUCT / n-ONLY ELIMINATION SEARCH")
print("=" * 78)

print()
print("random seed =", SEED)
print("R values    =", R_VALUES)
print("targets     =", len(targets))

print()
print("=" * 78)
print("REPRESENTATIVE TARGET")
print("=" * 78)

T0 = targets[len(targets) // 2]

print("p       =", T0["p"])
print("q       =", T0["q"])
print("n       =", T0["n"])
print("s=p+q   =", T0["s"])
print("n bits  =", T0["bits"])


# ============================================================================
# 1. FACTOR-SIDE CYCLOTOMIC STRUCTURE
# ============================================================================

print()
print("=" * 78)
print("1. FACTOR-SIDE CYCLOTOMIC STRUCTURE")
print("=" * 78)

p = T0["p"]
q = T0["q"]
n = T0["n"]
s = T0["s"]

Fp = F(p)
Fq = F(q)

print("F(p) =", Fp)
print("F(q) =", Fq)

print()
print("F(p)+F(q) =", Fp + Fq)
print("F(p)*F(q) =", Fp * Fq)
print("F(p)-F(q) =", Fp - Fq)

print()
print("Direct symmetric formulas:")

# F(p)+F(q)
sym_sum = s * s - 2 * n - s + 2

# F(p)F(q)
# (p²-p+1)(q²-q+1)
# = p²q² - p²q - pq² + p² + q² + pq - p - q + 1
# = n² - ns + (s²-2n) + n - s + 1
sym_prod = (
    n * n
    - n * s
    + s * s
    - n
    - s
    + 1
)

print("F(p)+F(q) symmetric formula =", sym_sum)
print("identity =", sym_sum == Fp + Fq)

print("F(p)F(q) symmetric formula =", sym_prod)
print("identity =", sym_prod == Fp * Fq)


# ============================================================================
# 2. ELIMINATING s FROM SIMPLE COMBINATIONS
# ============================================================================

print()
print("=" * 78)
print("2. SEARCH FOR s-ELIMINATING COMBINATIONS")
print("=" * 78)

print("""
We know:

    s = p+q
    n = pq

and

    (p-q)^2 = s^2 - 4n.

Any genuinely n-only factor invariant must eliminate s.

We search low-degree combinations of:

    s
    s^2
    n
    ns
    n^2
    F(p)+F(q)
    F(p)F(q)
    F(p)-F(q)
    (F(p)-F(q))^2
""")


# Build symbolic coefficient vectors in variables
# [1, s, s², n, ns, n²]
def factor_expression_vector(t):
    p = t["p"]
    q = t["q"]
    n = t["n"]
    s = t["s"]

    fp = F(p)
    fq = F(q)

    vals = {
        "1": 1,
        "s": s,
        "s2": s * s,
        "n": n,
        "ns": n * s,
        "n2": n * n,
        "Fp+Fq": fp + fq,
        "FpFq": fp * fq,
        "Fp-Fq": fp - fq,
        "diff2": (fp - fq) ** 2,
    }

    return vals


# Exact known eliminations.
print()
print("Known s-dependent identities:")

for name, expr in [
    ("F(p)+F(q)", Fp + Fq),
    ("F(p)F(q)", Fp * Fq),
    ("(F(p)-F(q))^2", (Fp - Fq) ** 2),
]:
    print(f"{name:25s} =", expr)


# ============================================================================
# 3. DISCRIMINANT SEARCH
# ============================================================================

print()
print("=" * 78)
print("3. DISCRIMINANT / FACTOR-DIFFERENCE SEARCH")
print("=" * 78)

print("""
The factor discriminant is

    Delta = s² - 4n = (p-q)².

We test whether cyclotomic factor expressions reduce to
simple functions of Delta and n.
""")

delta = s * s - 4 * n

print("Delta =", delta)
print("(p-q)^2 =", (p - q) ** 2)
print("identity =", delta == (p - q) ** 2)

print()
print("Useful exact identities:")

print("F(p)-F(q) =", (p - q) * (p + q - 1))
print("actual     =", Fp - Fq)
print("identity   =", (p - q) * (s - 1) == Fp - Fq)

print()
print("(F(p)-F(q))² = Delta*(s-1)²")

lhs = (Fp - Fq) ** 2
rhs = delta * (s - 1) ** 2

print("identity =", lhs == rhs)


# ============================================================================
# 4. MODULAR FACTOR-ROOT TEST
# ============================================================================

print()
print("=" * 78)
print("4. MODULAR FACTOR-ROOT TEST")
print("=" * 78)

print("""
For each auxiliary modulus m=F(r), test:

    F(p) mod m
    F(q) mod m
    p^6 mod m
    q^6 mod m

The key question is whether p or q themselves show a systematic
sixth-root relation modulo F(r).
""")

for r in R_VALUES[:9]:
    m = F(r)

    fp_mod = F(p) % m
    fq_mod = F(q) % m
    p6 = pow(p, 6, m)
    q6 = pow(q, 6, m)

    print(
        f"r={r:2d} m={m:6d} "
        f"F(p)={fp_mod:6d} "
        f"F(q)={fq_mod:6d} "
        f"p^6={p6:6d} "
        f"q^6={q6:6d}"
    )


# ============================================================================
# 5. PRODUCT OF FACTORS MOD F(r)
# ============================================================================

print()
print("=" * 78)
print("5. F(p)F(q) MOD F(r) VS n-ONLY EXPRESSIONS")
print("=" * 78)

print("""
We compare the factor-side quantity

    F(p)F(q) mod F(r)

against expressions computable from n alone:

    n
    n-1
    n+1
    n²
    n²-1
    n²-n+1
    n²+n+1

A match that survives across many independent targets is interesting.
""")

N_ONLY = {
    "n": lambda n: n,
    "n-1": lambda n: n - 1,
    "n+1": lambda n: n + 1,
    "n²": lambda n: n * n,
    "n²-1": lambda n: n * n - 1,
    "n²-n+1": lambda n: n * n - n + 1,
    "n²+n+1": lambda n: n * n + n + 1,
}

match_counts = {name: 0 for name in N_ONLY}
total_tests = 0

for t in targets:
    n = t["n"]
    p = t["p"]
    q = t["q"]

    fpq = F(p) * F(q)

    for r in R_VALUES:
        m = F(r)
        target_residue = fpq % m
        total_tests += 1

        for name, fn in N_ONLY.items():
            if target_residue == fn(n) % m:
                match_counts[name] += 1

print()

for name, count in match_counts.items():
    print(
        f"{name:12s} hits={count:5d}/{total_tests} "
        f"fraction={count/total_tests:.4f}"
    )


# ============================================================================
# 6. DIFFERENCE-SQUARE MODULAR SEARCH
# ============================================================================

print()
print("=" * 78)
print("6. (p-q)^2 / DELTA N-ONLY SEARCH")
print("=" * 78)

print("""
The obvious factor-sensitive quantity is

    Delta=(p-q)^2=s²-4n.

We test whether Delta modulo F(r) can accidentally coincide with
simple n-only expressions often enough to suggest a structural identity.
""")

DELTA_MATCHES = {
    name: 0
    for name in N_ONLY
}

for t in targets:
    p = t["p"]
    q = t["q"]
    n = t["n"]

    delta = (p - q) ** 2

    for r in R_VALUES:
        m = F(r)
        d = delta % m

        for name, fn in N_ONLY.items():
            if d == fn(n) % m:
                DELTA_MATCHES[name] += 1

for name, count in DELTA_MATCHES.items():
    print(
        f"{name:12s} hits={count:5d}/{total_tests} "
        f"fraction={count/total_tests:.4f}"
    )


# ============================================================================
# 7. SEARCH FOR SIMPLE UNIVERSAL FACTOR IDENTITIES
# ============================================================================

print()
print("=" * 78)
print("7. SIMPLE UNIVERSAL FACTOR IDENTITY SEARCH")
print("=" * 78)

print("""
We test whether expressions built from

    F(p), F(q), n

produce a universal constant modulo F(r).

Candidates:

    Fp+Fq-n
    FpFq-n²
    FpFq-n
    (Fp-Fq)²
    FpFq-(n²-1)
    FpFq-(n²-n+1)
    FpFq-(n²+n+1)
""")

candidate_functions = {
    "Fp+Fq-n":
        lambda fp, fq, n: fp + fq - n,

    "FpFq-n²":
        lambda fp, fq, n: fp * fq - n * n,

    "FpFq-n":
        lambda fp, fq, n: fp * fq - n,

    "(Fp-Fq)²":
        lambda fp, fq, n: (fp - fq) ** 2,

    "FpFq-(n²-1)":
        lambda fp, fq, n: fp * fq - (n * n - 1),

    "FpFq-(n²-n+1)":
        lambda fp, fq, n: fp * fq - (n * n - n + 1),

    "FpFq-(n²+n+1)":
        lambda fp, fq, n: fp * fq - (n * n + n + 1),
}

for name, fn in candidate_functions.items():
    values = []

    for t in targets:
        fp = F(t["p"])
        fq = F(t["q"])
        n = t["n"]

        values.append(fn(fp, fq, n))

    distinct = len(set(values))

    print(
        f"{name:28s} distinct={distinct:3d}"
    )


# ============================================================================
# 8. MODULAR COLLISION SEARCH
# ============================================================================

print()
print("=" * 78)
print("8. CROSS-TARGET MODULAR COLLISION SEARCH")
print("=" * 78)

print("""
For each modulus F(r), we group targets by n mod F(r).

If targets with the same n residue consistently have the same
factor-side residue

    F(p)F(q) mod F(r),

that would indicate that the factor expression is actually
determined by n modulo F(r).

We measure:

    ambiguity = number of different factor residues for one n residue.
""")

for r in R_VALUES[:9]:
    m = F(r)

    groups = {}

    for t in targets:
        nres = t["n"] % m
        fpqres = (F(t["p"]) * F(t["q"])) % m

        groups.setdefault(nres, set()).add(fpqres)

    ambiguous = {
        k: sorted(v)
        for k, v in groups.items()
        if len(v) > 1
    }

    print(
        f"r={r:2d} F={m:6d} "
        f"n_residues={len(groups):3d} "
        f"ambiguous={len(ambiguous):3d}"
    )

    if ambiguous:
        sample = list(ambiguous.items())[:3]
        print("  sample:", sample)


# ============================================================================
# 9. STRONGER COLLISION TEST: n + ALL AUXILIARY MODULI
# ============================================================================

print()
print("=" * 78)
print("9. MULTI-AUXILIARY n-RESIDUE COLLISION TEST")
print("=" * 78)

print("""
We now use the entire vector

    (n mod F(r1), ..., n mod F(rk))

and ask whether it determines

    (F(p)F(q) mod F(r1), ..., F(p)F(q) mod F(rk)).

If yes, that would be a much stronger n-only signal.
""")

for k in [2, 3, 5, 7, 9]:
    rs = R_VALUES[:k]

    groups = {}
    ambiguous = 0

    for t in targets:
        n_signature = tuple(t["n"] % F(r) for r in rs)

        factor_signature = tuple(
            (F(t["p"]) * F(t["q"])) % F(r)
            for r in rs
        )

        if n_signature not in groups:
            groups[n_signature] = set()

        groups[n_signature].add(factor_signature)

    ambiguous = sum(
        1 for vals in groups.values()
        if len(vals) > 1
    )

    print(
        f"k={k:2d} "
        f"groups={len(groups):3d} "
        f"ambiguous_groups={ambiguous:3d}"
    )


# ============================================================================
# 10. DIRECT SYMMETRIC ELIMINATION
# ============================================================================

print()
print("=" * 78)
print("10. DIRECT SYMMETRIC ELIMINATION OF s")
print("=" * 78)

print("""
We explicitly solve the factor-side formulas in terms of

    n=pq
    s=p+q.

For F(x)=x²-x+1:

    S = F(p)+F(q)
      = s² - s - 2n + 2

    P = F(p)F(q)
      = n² - ns + s² - n - s + 1

Using S to eliminate s²:

    P-S = n² - ns + s - n - 1

    P-S = n²-n-1 - s(n-1)

Therefore, unless s disappears through an additional relation,
the product still contains p+q.

This is the exact obstruction we test numerically.
""")

for t in targets[:3]:
    p = t["p"]
    q = t["q"]
    n = t["n"]
    s = t["s"]

    S = F(p) + F(q)
    P = F(p) * F(q)

    reduced = n * n - n - 1 - s * (n - 1)

    print()
    print("target n =", n)
    print("P-S      =", P - S)
    print("reduced  =", reduced)
    print("identity =", P - S == reduced)


# ============================================================================
# 11. SEARCH FOR CANCELLATION OF s
# ============================================================================

print()
print("=" * 78)
print("11. LINEAR COMBINATION SEARCH FOR s CANCELLATION")
print("=" * 78)

print("""
Search expressions

    a*(F(p)+F(q))
  + b*(F(p)F(q))
  + c*n²
  + d*n
  + e

with small coefficients.

A genuine n-only identity exists if the coefficient of s and s²
vanishes symbolically.

We perform the symbolic coefficient test directly.
""")

# Symbolic representation in [1, s, s², n, ns, n²]
# S = s² - s - 2n + 2
# P = s² - ns - n - s + n² + 1

S_vec = {
    "1": 2,
    "s": -1,
    "s2": 1,
    "n": -2,
    "ns": 0,
    "n2": 0,
}

P_vec = {
    "1": 1,
    "s": -1,
    "s2": 1,
    "n": -1,
    "ns": -1,
    "n2": 1,
}

found = []

for a in range(-3, 4):
    for b in range(-3, 4):
        for c in range(-3, 4):
            for d in range(-3, 4):
                for e in range(-3, 4):

                    if (a, b, c, d, e) == (0, 0, 0, 0, 0):
                        continue

                    coeff = {
                        key:
                        a * S_vec.get(key, 0)
                        + b * P_vec.get(key, 0)
                        + c * (1 if key == "n2" else 0)
                        + d * (1 if key == "n" else 0)
                        + e * (1 if key == "1" else 0)
                        for key in ["1", "s", "s2", "n", "ns", "n2"]
                    }

                    # Require elimination of all s-containing terms.
                    if coeff["s"] == 0 and coeff["s2"] == 0 and coeff["ns"] == 0:
                        found.append((a, b, c, d, e, coeff))

print("s-free linear combinations found =", len(found))

for item in found[:PRINT_LIMIT]:
    a, b, c, d, e, coeff = item
    print(
        f"({a})S + ({b})P + ({c})n² + ({d})n + ({e})"
    )


# ============================================================================
# 12. MODULAR EXISTENCE TEST
# ============================================================================

print()
print("=" * 78)
print("12. MODULAR EXISTENCE OF FACTOR ROOT PAIRS")
print("=" * 78)

print("""
For each m=F(r), define

    Z_m = {x mod m : F(x)=0 mod m}.

If p and q were constrained by the cyclotomic auxiliary structure,
then their residues might belong to Z_m.

We measure how often the actual p,q satisfy this condition.
""")

for r in R_VALUES[:9]:
    m = F(r)

    roots = [x for x in range(m) if F(x) % m == 0]

    p_root = p % m in roots
    q_root = q % m in roots

    print(
        f"r={r:2d} m={m:6d} "
        f"roots={roots[:12]} "
        f"p_root={p_root} q_root={q_root}"
    )


# ============================================================================
# 13. n-ONLY PRODUCT SET
# ============================================================================

print()
print("=" * 78)
print("13. CAN n MOD F(r) DETERMINE A FACTOR-PRODUCT CLASS?")
print("=" * 78)

print("""
For each m=F(r), enumerate products

    x*y mod m

where x,y range over all residue classes satisfying the
sixth-root relation

    x^6 = 1 mod m
    y^6 = 1 mod m.

Then compare n mod m against this product set.

This is a purely n,r computation.
No p,q are used in constructing the set.
""")

for r in R_VALUES[:9]:
    m = F(r)

    roots6 = [
        x for x in range(m)
        if pow(x, 6, m) == 1
    ]

    product_set = {
        (x * y) % m
        for x in roots6
        for y in roots6
    }

    actual_n = n % m

    print(
        f"r={r:2d} m={m:6d} "
        f"sixth_roots={len(roots6):3d} "
        f"product_classes={len(product_set):4d} "
        f"n_mod_m={actual_n:6d} "
        f"in_product_set={actual_n in product_set}"
    )


# ============================================================================
# 14. MULTI-AUXILIARY PRODUCT CLASS COLLAPSE
# ============================================================================

print()
print("=" * 78)
print("14. MULTI-AUXILIARY PRODUCT-CLASS COLLAPSE")
print("=" * 78)

print("""
We combine the sixth-root product constraints over several F(r).

The hope would be:

    n alone
      ->
    a tiny number of possible factor residue pairs.

We measure the number of possible product signatures.
""")

for k in [1, 2, 3, 4, 5]:

    rs = R_VALUES[:k]

    possible = [set()]

    for r in rs:
        m = F(r)

        roots6 = [
            x for x in range(m)
            if pow(x, 6, m) == 1
        ]

        products = {
            (x * y) % m
            for x in roots6
            for y in roots6
        }

        possible.append(products)

    signature_count = 1

    for r in rs:
        m = F(r)
        roots6 = [
            x for x in range(m)
            if pow(x, 6, m) == 1
        ]

        products = {
            (x * y) % m
            for x in roots6
            for y in roots6
        }

        signature_count *= len(products)

    print(
        f"k={k:2d} "
        f"moduli={[F(r) for r in rs]} "
        f"independent_product_classes={signature_count}"
    )


# ============================================================================
# 15. FINAL CLASSIFICATION
# ============================================================================

print()
print("=" * 78)
print("15. FINAL CLASSIFICATION")
print("=" * 78)

print("""
A. FACTOR-SIDE STRUCTURE
------------------------
F(p)+F(q) and F(p)F(q) are symmetric in p,q.

But both generally depend on

    s=p+q

in addition to n=pq.

B. DISCRIMINANT
---------------
(p-q)^2 = s²-4n

also retains the unknown s.

C. CYCLOTOMIC CONDITION
-----------------------
F(r)=Phi_6(r) is intrinsically associated with sixth roots of unity.

This does not by itself impose

    F(p)=0 mod F(r)

or

    F(q)=0 mod F(r).

D. n-ONLY TEST
--------------
We test whether n mod F(r) determines factor-side cyclotomic
products.

The important measurement is the collision/ambiguity count.

E. BREAKTHROUGH CONDITION
-------------------------
A result would be interesting only if:

    n mod F(r)
       ->
    unique/small factor-side class
       ->
    additional constraint on u or A.

If the mapping remains highly ambiguous, the cyclotomic
factor route is not sufficient.

F. MOST IMPORTANT ALGEBRAIC OBSTRUCTION
---------------------------------------
The exact identities are:

    F(p)+F(q) = s²-s-2n+2

and

    F(p)F(q)
      = n²-n+1 - s(n+1) + s².

Thus the factor-side information carries s=p+q.

The next question is whether some combination of several
cyclotomic evaluations cancels s completely while leaving
a nontrivial n-dependent residue.

==============================================================================
SUMMARY
==============================================================================

Targets                    = {len(targets)}
Auxiliary r values         = {R_VALUES}

This experiment deliberately moves away from generic polynomial
search and attacks the p*q structure directly.

No CSV files are produced.
Only printed diagnostics are generated.

==============================================================================
DONE
==============================================================================
""".format(
    len(targets),
    R_VALUES
))

