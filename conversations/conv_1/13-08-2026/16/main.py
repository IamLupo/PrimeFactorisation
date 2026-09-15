#!/usr/bin/env python3

"""
==============================================================================
KAPPA NEXT EXPERIMENT
CYCLOTOMIC CHARACTER / FACTOR-CLASS SIGNATURE SEARCH
==============================================================================

Goal
----

The previous adversarial experiment established a strong negative result:

    n mod F(r)  does NOT determine  s mod F(r)

because unit transformations

    (x,y) -> (x*t, y*t^(-1))

preserve xy but generally change x+y.

So we stop trying to recover s directly.

This experiment instead asks:

    Does the collection of cyclotomic / cubic residue information
    attached to the factor residues contain an n-only invariant?

Core identity:

    F(x) = x^2 - x + 1

and

    x^3 + 1 = (x+1) F(x).

Therefore, modulo m = F(r), the condition

    F(x) = 0 mod m

is equivalent, for x != -1, to

    x^3 = -1 mod m.

We investigate:

    x^3 mod m
    (x+1)^3 mod m
    F(x) mod m
    Legendre/Jacobi-style character information where defined
    discrete cubic-class signatures
    factor-pair signatures

The critical control is adversarial:

For the SAME product residue n mod m, construct many unit factor pairs

    x*y = n mod m

and measure whether any proposed signature is invariant across
all such pairs.

If a signature varies freely under the unit action, it is useless
as an n-only invariant.

If some signature survives all adversarial transformations, it
becomes a candidate structural invariant.

NO CSV.
ONLY PRINTS.
==============================================================================
"""

from __future__ import annotations

import math
import random
from collections import defaultdict, Counter
from itertools import combinations


# ============================================================================
# CONFIG
# ============================================================================

SEED = 20260814

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47
]

TARGETS = 24

PRIME_BITS = 22

PAIR_SAMPLES = 250

COMBINED_K = [2, 3, 4, 5]

random.seed(SEED)


# ============================================================================
# BASIC FUNCTIONS
# ============================================================================

def F(x: int) -> int:
    return x * x - x + 1


def is_prime(n: int) -> bool:
    if n < 2:
        return False

    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    # More than enough for our ~22-bit primes.
    for a in [2, 3, 5, 7, 11]:
        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                break
        else:
            return False

    return True


def random_prime(bits: int) -> int:
    while True:
        x = random.getrandbits(bits)
        x |= (1 << (bits - 1))
        x |= 1

        if is_prime(x):
            return x


def unit_residues(m: int):
    return [x for x in range(m) if math.gcd(x, m) == 1]


def sixth_roots(m: int):
    return [x for x in unit_residues(m) if pow(x, 6, m) == 1]


def cubic_minus_one_roots(m: int):
    return [x for x in range(m) if pow(x, 3, m) == (-1) % m]


def gcd_safe(a: int, m: int) -> int:
    return math.gcd(a % m, m)


# ============================================================================
# JACOBI SYMBOL
# ============================================================================

def jacobi(a: int, n: int) -> int:
    """
    Exact Jacobi symbol (a/n), n odd positive.
    Returns -1, 0, +1.
    """

    if n <= 0 or n % 2 == 0:
        raise ValueError("Jacobi denominator must be positive odd")

    a %= n
    result = 1

    while a:
        while a % 2 == 0:
            a //= 2
            r = n % 8

            if r == 3 or r == 5:
                result = -result

        a, n = n, a

        if a % 4 == 3 and n % 4 == 3:
            result = -result

        a %= n

    return result if n == 1 else 0


# ============================================================================
# CYCLOTOMIC SIGNATURE
# ============================================================================

def cyclotomic_signature(x: int, m: int):
    """
    A deliberately redundant signature.

    We want to see which components survive adversarial
    product-preserving transformations.
    """

    x %= m

    xp1 = (x + 1) % m

    return (
        x,
        pow(x, 2, m),
        pow(x, 3, m),
        pow(x, 6, m),
        F(x) % m,
        xp1,
        pow(xp1, 2, m),
        pow(xp1, 3, m),
        gcd_safe(x + 1, m),
    )


def pair_signature(x: int, y: int, m: int):
    """
    Symmetric factor-pair signature.
    """

    x %= m
    y %= m

    fx = F(x) % m
    fy = F(y) % m

    return {
        "Fsum": (fx + fy) % m,
        "Fprod": (fx * fy) % m,
        "Fdiff2": (fx - fy) ** 2 % m,

        "cube_sum": (pow(x, 3, m) + pow(y, 3, m)) % m,
        "cube_prod": (pow(x, 3, m) * pow(y, 3, m)) % m,

        "xp1_yp1": ((x + 1) * (y + 1)) % m,

        "xp1_sum": (x + 1 + y + 1) % m,

        "gcd_xp1": gcd_safe(x + 1, m),
        "gcd_yp1": gcd_safe(y + 1, m),
    }


# ============================================================================
# PRODUCT-PRESERVING ORBIT
# ============================================================================

def orbit_pairs(n_mod: int, m: int):
    """
    Generate all unit factor pairs

        x*y = n_mod mod m

    whenever n_mod is itself a unit.

    For every unit t:

        x=t
        y=n*t^{-1}

    This is the exact adversarial symmetry.
    """

    if math.gcd(n_mod, m) != 1:
        return []

    inv_n = n_mod % m

    result = []

    for x in unit_residues(m):
        y = (n_mod * pow(x, -1, m)) % m

        result.append((x, y))

    return result


# ============================================================================
# REPRESENTATIVE TARGETS
# ============================================================================

def generate_targets(count: int):
    targets = []

    while len(targets) < count:
        p = random_prime(PRIME_BITS)
        q = random_prime(PRIME_BITS)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q
        s = p + q
        delta = (p - q) ** 2

        targets.append({
            "p": p,
            "q": q,
            "n": n,
            "s": s,
            "delta": delta,
        })

    return targets


# ============================================================================
# SECTION 1
# ============================================================================

def section_header(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def print_target(t):
    print(f"p       = {t['p']}")
    print(f"q       = {t['q']}")
    print(f"n       = {t['n']}")
    print(f"s       = {t['s']}")
    print(f"Delta   = {t['delta']}")
    print(f"n bits  = {t['n'].bit_length()}")


# ============================================================================
# SECTION 2: ROOT STRUCTURE
# ============================================================================

def root_structure():
    section_header("1. CYCLOTOMIC ROOT STRUCTURE")

    for r in R_VALUES:
        m = F(r)

        roots_f = [
            x for x in range(m)
            if F(x) % m == 0
        ]

        roots_cubic = cubic_minus_one_roots(m)

        print(
            f"r={r:2d} "
            f"m={m:6d} "
            f"F-roots={len(roots_f):4d} "
            f"x^3=-1 roots={len(roots_cubic):4d}"
        )

        if m <= 100:
            print(f"  F-roots       = {roots_f}")
            print(f"  cubic roots   = {roots_cubic}")

        overlap = set(roots_f) & set(roots_cubic)

        print(
            f"  exact overlap = {len(overlap):4d} "
            f"identity_check={all(pow(x, 3, m) == (-1) % m for x in roots_f)}"
        )


# ============================================================================
# SECTION 3: ACTUAL FACTOR SIGNATURES
# ============================================================================

def actual_factor_signatures(targets):
    section_header("2. ACTUAL FACTOR CYCLOTOMIC SIGNATURES")

    for idx, t in enumerate(targets[:8], 1):
        print()
        print(f"TARGET {idx}")

        for r in R_VALUES[:9]:
            m = F(r)

            p = t["p"] % m
            q = t["q"] % m

            sig_p = cyclotomic_signature(p, m)
            sig_q = cyclotomic_signature(q, m)

            print(
                f"r={r:2d} m={m:5d} "
                f"p={p:5d} q={q:5d} "
                f"Fp={F(p)%m:5d} Fq={F(q)%m:5d} "
                f"p3={pow(p,3,m):5d} "
                f"q3={pow(q,3,m):5d}"
            )


# ============================================================================
# SECTION 4: ADVERSARIAL INVARIANT TEST
# ============================================================================

def adversarial_invariant_test():
    section_header("3. ADVERSARIAL PRODUCT-ORBIT INVARIANT TEST")

    """
    For each modulus choose several product residues n.

    For every pair (x,y) with xy=n, compare candidate signatures.

    A candidate is useful only if it remains constant throughout
    the complete product-preserving orbit.
    """

    candidate_names = [
        "Fsum",
        "Fprod",
        "Fdiff2",
        "cube_sum",
        "cube_prod",
        "xp1_yp1",
        "xp1_sum",
        "gcd_xp1",
        "gcd_yp1",
    ]

    global_stats = {
        name: {"tested": 0, "invariant": 0}
        for name in candidate_names
    }

    for r in R_VALUES:
        m = F(r)

        units = unit_residues(m)

        if not units:
            continue

        sample_products = []

        for _ in range(min(PAIR_SAMPLES, len(units))):
            x = random.choice(units)
            y = random.choice(units)
            sample_products.append((x * y) % m)

        sample_products = sorted(set(sample_products))

        print()
        print(f"r={r:2d} m={m:6d} product_classes_tested={len(sample_products)}")

        for n_mod in sample_products:

            pairs = orbit_pairs(n_mod, m)

            if not pairs:
                continue

            if len(pairs) > PAIR_SAMPLES:
                pairs = random.sample(pairs, PAIR_SAMPLES)

            values = {
                name: set()
                for name in candidate_names
            }

            for x, y in pairs:
                ps = pair_signature(x, y, m)

                for name in candidate_names:
                    values[name].add(ps[name])

            for name in candidate_names:
                global_stats[name]["tested"] += 1

                if len(values[name]) == 1:
                    global_stats[name]["invariant"] += 1

        print("  orbit-invariant fraction:")

        for name in candidate_names:
            st = global_stats[name]

            frac = (
                st["invariant"] / st["tested"]
                if st["tested"]
                else 0.0
            )

            print(
                f"    {name:12s} "
                f"{st['invariant']:5d}/{st['tested']:5d} "
                f"{frac:8.4f}"
            )


# ============================================================================
# SECTION 5: CHARACTER SEARCH
# ============================================================================

def character_search():
    section_header("4. JACOBI / QUADRATIC CHARACTER SEARCH")

    """
    We test whether quadratic-character information survives the
    product-preserving action.

    For odd m we factor m into prime powers and compute Jacobi
    symbols where meaningful.

    The important point is that

        chi(x) * chi(y) = chi(xy)

    for multiplicative characters.

    Hence any multiplicative character of a factor is automatically
    an n-only invariant after multiplying the two factor characters.

    We explicitly verify this and search for less trivial
    pair expressions.
    """

    for r in R_VALUES:
        m = F(r)

        if m % 2 == 0:
            continue

        units = unit_residues(m)

        if not units:
            continue

        good = 0
        total = 0

        for _ in range(min(500, len(units) ** 2)):
            x = random.choice(units)
            y = random.choice(units)

            n = (x * y) % m

            chi_x = jacobi(x, m)
            chi_y = jacobi(y, m)
            chi_n = jacobi(n, m)

            total += 1

            if chi_x * chi_y == chi_n:
                good += 1

        print(
            f"r={r:2d} m={m:6d} "
            f"Jacobi multiplicativity={good}/{total} "
            f"fraction={good/total if total else 0:.4f}"
        )


# ============================================================================
# SECTION 6: CUBIC CHARACTER PROXY
# ============================================================================

def cubic_class(x: int, m: int):
    """
    Return the multiplicative order of x modulo m when x is a unit.

    We use this as a robust finite-group class proxy.

    For a unit modulo m, the order divides phi(m).
    """

    if math.gcd(x, m) != 1:
        return 0

    phi = 0
    # Exact phi by direct counting is acceptable because m <= 2163.
    for z in range(1, m + 1):
        if math.gcd(z, m) == 1:
            phi += 1

    order = phi

    for p in range(2, int(math.isqrt(order)) + 1):
        if order % p == 0:
            while order % p == 0 and pow(x, order // p, m) == 1:
                order //= p

    return order


def cubic_structure_search():
    section_header("5. MULTIPLICATIVE ORDER / CUBIC-CLASS SEARCH")

    for r in R_VALUES:
        m = F(r)

        units = unit_residues(m)

        if not units:
            continue

        order_counts = Counter()

        for x in units:
            order_counts[cubic_class(x, m)] += 1

        print(
            f"r={r:2d} m={m:6d} "
            f"unit_count={len(units):5d} "
            f"orders={dict(sorted(order_counts.items()))}"
        )


# ============================================================================
# SECTION 7: FACTOR-PAIR ORDER PRODUCT
# ============================================================================

def order_pair_test():
    section_header("6. FACTOR-PAIR ORDER RELATION")

    """
    If x*y=n, multiplicative order data satisfies constraints,
    but order(x) and order(y) individually need not be determined
    by n.

    We measure exactly how ambiguous these values are for fixed n.
    """

    for r in R_VALUES[:12]:
        m = F(r)

        units = unit_residues(m)

        if len(units) < 2:
            continue

        n_mod = random.choice(units)

        pairs = orbit_pairs(n_mod, m)

        if len(pairs) > PAIR_SAMPLES:
            pairs = random.sample(pairs, PAIR_SAMPLES)

        orders = set()
        gcd_orders = set()
        lcm_orders = set()

        for x, y in pairs:
            ox = cubic_class(x, m)
            oy = cubic_class(y, m)

            orders.add((ox, oy))
            gcd_orders.add(math.gcd(ox, oy))
            lcm_orders.add(math.lcm(ox, oy))

        print(
            f"r={r:2d} m={m:6d} "
            f"n={n_mod:6d} "
            f"order_pairs={len(orders):5d} "
            f"gcd_orders={sorted(gcd_orders)} "
            f"lcm_orders={sorted(lcm_orders)}"
        )


# ============================================================================
# SECTION 8: SPECIAL CYCLOTOMIC EVENT
# ============================================================================

def special_event_search(targets):
    section_header("7. SPECIAL EVENT SEARCH: F(p)=0 OR F(q)=0")

    """
    Test whether actual factor residues frequently enter the
    cyclotomic root classes modulo auxiliary F(r).

    This is deliberately compared against the exact expected
    root density among units.

    If the actual factors show a strong excess, that would be
    a meaningful signal.
    """

    for r in R_VALUES:
        m = F(r)

        roots = {
            x for x in unit_residues(m)
            if F(x) % m == 0
        }

        hits_p = 0
        hits_q = 0
        total = len(targets) * 2

        for t in targets:
            if t["p"] % m in roots:
                hits_p += 1

            if t["q"] % m in roots:
                hits_q += 1

        root_density = (
            len(roots) / len(unit_residues(m))
            if unit_residues(m)
            else 0
        )

        print(
            f"r={r:2d} m={m:6d} "
            f"roots={len(roots):3d} "
            f"unit_density={root_density:.6f} "
            f"p_hits={hits_p:3d}/{len(targets)} "
            f"q_hits={hits_q:3d}/{len(targets)}"
        )


# ============================================================================
# SECTION 9: COMBINED CHARACTER SIGNATURE
# ============================================================================

def combined_character_signature(x, moduli):
    """
    n-only candidate:

        product over moduli of multiplicative character data.

    We intentionally use only information that is multiplicative
    in x, so that pair products may collapse to n-only values.
    """

    sig = []

    for m in moduli:
        if math.gcd(x, m) != 1:
            sig.append(None)
            continue

        sig.append((
            jacobi(x, m) if m % 2 else None,
            cubic_class(x, m),
            pow(x, 3, m),
            pow(x, 6, m),
        ))

    return tuple(sig)


def combined_signature_test():
    section_header("8. MULTI-MODULUS CHARACTER COLLISION TEST")

    """
    This is the main new direction.

    For a factor pair (x,y), compute:

        character(x) * character(y)

    and compare with the corresponding signature of

        n = x*y.

    Anything that agrees automatically is an n-only invariant.

    We then search for pair signatures that are NOT obviously
    multiplicative and see whether they nevertheless collapse.
    """

    for k in COMBINED_K:
        moduli = [F(r) for r in R_VALUES[:k]]

        print()
        print(
            f"k={k} "
            f"moduli={moduli}"
        )

        total = 0
        multiplicative_ok = 0

        nontrivial_classes = defaultdict(set)

        for _ in range(300):

            x = random.randrange(2, 100000)
            y = random.randrange(2, 100000)

            if any(math.gcd(x, m) != 1 or math.gcd(y, m) != 1
                   for m in moduli):
                continue

            n = x * y

            sx = combined_character_signature(x, moduli)
            sy = combined_character_signature(y, moduli)
            sn = combined_character_signature(n, moduli)

            total += 1

            # Only test the multiplicative-order component.
            order_product = []

            valid = True

            for a, b, c in zip(sx, sy, sn):

                if a[1] == 0 or b[1] == 0 or c[1] == 0:
                    valid = False
                    break

                order_product.append(
                    (a[1], b[1], c[1])
                )

            if valid:
                # n-only consistency check:
                # order(n) must divide lcm(order(x), order(y)).
                ok = True

                for ox, oy, on in order_product:
                    if math.lcm(ox, oy) % on != 0:
                        ok = False

                if ok:
                    multiplicative_ok += 1

            nontrivial_classes[
                tuple(z[1] for z in sx)
            ].add(tuple(z[1] for z in sy))

        print(
            f"  samples={total}"
        )

        print(
            f"  order divisibility consistency="
            f"{multiplicative_ok}/{total}"
        )

        print(
            f"  x-order classes="
            f"{len(nontrivial_classes)}"
        )


# ============================================================================
# SECTION 10: SEARCH FOR AN ACTUAL FACTORIZATION SIGNAL
# ============================================================================

def factor_signal_search(targets):
    section_header("9. ACTUAL TARGET N-ONLY SIGNAL SEARCH")

    """
    For each actual target, compute n-only signatures.

    We deliberately DO NOT use p,q,s,Delta to construct the
    candidate signature.

    We ask whether the signature separates targets unusually well.

    This is not claiming that separation factors n.
    It is simply a screening test for useful arithmetic structure.
    """

    for k in [1, 2, 3, 4, 5]:

        moduli = [F(r) for r in R_VALUES[:k]]

        signatures = defaultdict(list)

        for i, t in enumerate(targets):

            sig = []

            for m in moduli:
                nm = t["n"] % m

                sig.append((
                    nm,
                    pow(nm, 2, m),
                    pow(nm, 3, m),
                    jacobi(nm, m) if m % 2 else None,
                ))

            signatures[tuple(sig)].append(i)

        ambiguous = [
            v for v in signatures.values()
            if len(v) > 1
        ]

        print(
            f"k={k} "
            f"moduli={moduli} "
            f"groups={len(signatures):3d} "
            f"ambiguous_groups={len(ambiguous):3d} "
            f"max_group={max(map(len, signatures.values()))}"
        )


# ============================================================================
# SECTION 11: THE KEY CONTROL
# ============================================================================

def key_control():
    section_header("10. KEY ALGEBRAIC CONTROL")

    print(
        """
For every unit t modulo M:

    (x,y) -> (xt, y*t^(-1))

preserves:

    xy mod M.

Therefore any proposed factor-side quantity Q(x,y) is a valid
n-only invariant only if

    Q(xt, y*t^(-1)) = Q(x,y)

for every admissible t.

This experiment therefore distinguishes two cases:

    CASE A
    ------
    Q is invariant under the complete unit action.

    Then Q may genuinely be a function of n.

    CASE B
    ------
    Q changes under the unit action.

    Then Q cannot be an n-only invariant.

The previous experiment showed:

    s
    F(x)+F(y)
    F(x)-F(y)
    Delta

all fail this adversarial test.

The present experiment asks whether multiplicative cyclotomic
character data survives where additive symmetric data failed.
"""
    )


# ============================================================================
# SECTION 12: FINAL CLASSIFICATION
# ============================================================================

def final_classification():
    section_header("11. FINAL CLASSIFICATION")

    print(
        """
A. IF MULTIPLICATIVE CHARACTER DATA COLLAPSES TO n
---------------------------------------------------

This is expected for genuine multiplicative characters.

For example:

    chi(x) chi(y) = chi(xy) = chi(n).

That would be an n-only invariant, but not automatically a
factorization mechanism.

B. MORE INTERESTING CASE
------------------------

A genuinely interesting result would be a factor-side quantity
which:

    1. is not trivially multiplicative,
    2. survives the complete unit-action control,
    3. persists across many F(r),
    4. distinguishes actual factor classes,
    5. gives additional constraints beyond n mod M.

C. BREAKTHROUGH CONDITION
-------------------------

The next useful threshold is NOT merely:

    "we found another invariant."

It is:

    n
      ->
    restricted factor residue classes
      ->
    restricted s or Delta
      ->
    small number of integer candidates
      ->
    exact factorization.

D. IMPORTANT NEGATIVE OUTCOME
-----------------------------

If every nontrivial candidate varies under

    (x,y) -> (xt,y*t^(-1)),

then the auxiliary cyclotomic moduli are behaving exactly like
generic residue moduli.

At that point the next experiment should move away from residue
classes entirely and investigate a different algebraic object,
such as resultants, norms, or higher cyclotomic polynomials.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    section_header(
        "KAPPA NEXT EXPERIMENT\n"
        "CYCLOTOMIC CHARACTER / FACTOR-CLASS SIGNATURE SEARCH"
    )

    print(f"random seed = {SEED}")
    print(f"R values    = {R_VALUES}")
    print(f"targets     = {TARGETS}")
    print("CSV output  = NONE")

    targets = generate_targets(TARGETS)

    section_header("REPRESENTATIVE TARGET")
    print_target(targets[0])

    root_structure()
    actual_factor_signatures(targets)
    adversarial_invariant_test()
    character_search()
    cubic_structure_search()
    order_pair_test()
    special_event_search(targets)
    combined_signature_test()
    factor_signal_search(targets)
    key_control()
    final_classification()

    section_header("DONE")


if __name__ == "__main__":
    main()

