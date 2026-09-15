#!/usr/bin/env python3

"""
==============================================================================
KAPPA NEXT EXPERIMENT
FULL ORBIT VS ACTUAL FACTOR DISTRIBUTION
ADVERSARIAL CYCLOTOMIC SIGNATURE TEST
==============================================================================

Goal
----
The previous experiment found that combining several F(r) moduli can
collapse the observed 24 targets into unique groups.

That is NOT yet evidence of a factorization mechanism.

This experiment asks the stronger question:

    Given n = p*q mod M,

    among ALL admissible factor pairs

        (x, y) with x*y == n (mod M),

    does the actual pair

        (p mod M, q mod M)

    occupy an unusually small / special subset?

We therefore exhaustively enumerate the complete unit orbit:

    (x, y) -> (x*t, y*t^(-1))

and compare the actual prime-factor pair against every admissible pair.

No CSV output.
Prints only.

==============================================================================
"""

from math import gcd
from collections import Counter, defaultdict
from itertools import combinations

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

SEED = 20260814

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17,
    19, 23, 29, 31, 37, 41, 43, 47
]

TARGET_COUNT = 24

# Use the representative target plus deterministically generated targets.
REPRESENTATIVE_P = 3622631
REPRESENTATIVE_Q = 3922519

# Number of printed orbit witnesses per modulus.
MAX_WITNESSES = 8

# Keep exhaustive orbit calculations for moduli up to this size.
# All configured moduli are small enough, but this guard protects future use.
MAX_EXHAUSTIVE_M = 10000

# Combined modulus stages.
COMBINED_STAGES = [
    [2],
    [2, 3],
    [2, 3, 5],
    [2, 3, 5, 7],
    [2, 3, 5, 7, 11],
]


# ---------------------------------------------------------------------------
# F(r)
# ---------------------------------------------------------------------------

def F(r):
    return r * r + r + 1


# ---------------------------------------------------------------------------
# Deterministic primality / prime generation
# ---------------------------------------------------------------------------

def is_prime(n):
    if n < 2:
        return False

    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)

    if n in small:
        return True

    for p in small:
        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    # Deterministic for our small primes.
    for a in (2, 3, 5, 7, 11, 13):
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


def next_prime(n):
    if n <= 2:
        return 2

    if n % 2 == 0:
        n += 1

    while not is_prime(n):
        n += 2

    return n


def generate_targets(count):
    """
    Deterministically generate distinct prime pairs.
    The representative target is TARGET 1.
    """
    targets = [(REPRESENTATIVE_P, REPRESENTATIVE_Q)]

    p = 8000003
    q = 9000001

    while len(targets) < count:
        p = next_prime(p)
        q = next_prime(q)

        if p == q:
            q = next_prime(q + 2)

        if (p, q) not in targets and (q, p) not in targets:
            targets.append((p, q))

        p += 7919
        q += 104729

    return targets


# ---------------------------------------------------------------------------
# Unit-group helpers
# ---------------------------------------------------------------------------

def units(m):
    return [x for x in range(1, m) if gcd(x, m) == 1]


def unit_inverse(x, m):
    return pow(x, -1, m)


def multiplicative_order(x, m):
    if gcd(x, m) != 1:
        return None

    y = 1

    for k in range(1, 100000):
        y = (y * x) % m

        if y == 1:
            return k

    raise RuntimeError("order search exceeded bound")


def jacobi_symbol(a, n):
    """
    General Jacobi symbol.
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


# ---------------------------------------------------------------------------
# Cyclotomic helpers
# ---------------------------------------------------------------------------

def cyclotomic_F(x):
    """
    x^2 + x + 1
    """
    return x * x + x + 1


def cube_class(x, m):
    return pow(x, 3, m)


def cubic_root_flag(x, m):
    return pow(x, 3, m) == (m - 1) % m


def F_zero_flag(x, m):
    return cyclotomic_F(x) % m == 0


def diff2(x, y, m):
    return ((x - y) ** 2) % m


def cube_sum(x, y, m):
    return (pow(x, 3, m) + pow(y, 3, m)) % m


def cube_diff(x, y, m):
    return (pow(x, 3, m) - pow(y, 3, m)) % m


def xp1_yp1(x, y, m):
    return ((x + 1) * (y + 1)) % m


def xp1_sum(x, y, m):
    return ((x + 1) * (x + y)) % m


def gcd_xp1(x, m):
    return gcd(x + 1, m)


def gcd_yp1(y, m):
    return gcd(y + 1, m)


# ---------------------------------------------------------------------------
# Signature definitions
# ---------------------------------------------------------------------------

def signatures(x, y, m):
    """
    Candidate factor-side signatures.

    Some are intentionally redundant. That is useful because we want to
    distinguish:

        purely multiplicative invariants
        symmetric but non-invariant quantities
        cyclotomic quantities
        order-based quantities
    """

    ox = multiplicative_order(x, m)
    oy = multiplicative_order(y, m)

    return {
        "sum": (x + y) % m,

        "diff2": diff2(x, y, m),

        "cube_sum": cube_sum(x, y, m),

        "cube_diff": cube_diff(x, y, m),

        "cube_prod": (pow(x, 3, m) * pow(y, 3, m)) % m,

        "Fsum": (
            cyclotomic_F(x) + cyclotomic_F(y)
        ) % m,

        "Fprod": (
            cyclotomic_F(x) * cyclotomic_F(y)
        ) % m,

        "Fdiff": (
            cyclotomic_F(x) - cyclotomic_F(y)
        ) % m,

        "xp1_yp1": xp1_yp1(x, y, m),

        "xp1_sum": xp1_sum(x, y, m),

        "gcd_xp1": gcd_xp1(x, m),

        "gcd_yp1": gcd_yp1(y, m),

        "order_pair": (
            min(ox, oy),
            max(ox, oy),
        ),

        "order_product": (
            ox * oy
        ),

        "order_gcd": gcd(ox, oy),

        "order_lcm": (ox * oy) // gcd(ox, oy),

        "cube_class_pair": (
            pow(x, 3, m),
            pow(y, 3, m),
        ),

        "F_zero_pair": (
            F_zero_flag(x, m),
            F_zero_flag(y, m),
        ),

        "Jacobi_pair": (
            jacobi_symbol(x, m),
            jacobi_symbol(y, m),
        ),
    }


# ---------------------------------------------------------------------------
# Exact product orbit
# ---------------------------------------------------------------------------

def product_orbit(n_mod, m):
    """
    Enumerate all ordered unit factor pairs:

        x*y == n_mod (mod m).

    Since y = n*x^{-1}, every unit x produces exactly one y.
    """

    U = units(m)
    result = []

    for x in U:
        y = (n_mod * unit_inverse(x, m)) % m

        if gcd(y, m) != 1:
            raise AssertionError("constructed y is not a unit")

        if (x * y) % m != n_mod:
            raise AssertionError("orbit product failure")

        result.append((x, y))

    return result


# ---------------------------------------------------------------------------
# Actual target representation
# ---------------------------------------------------------------------------

def target_residue_pair(p, q, m):
    return p % m, q % m


# ---------------------------------------------------------------------------
# Signature rarity
# ---------------------------------------------------------------------------

def signature_rarity(orbit, actual_pair, m, signature_name):
    """
    Returns:

        actual signature
        number of orbit pairs sharing it
        orbit size
        fraction
        rank-like percentile
    """

    actual_sig = signatures(
        actual_pair[0],
        actual_pair[1],
        m
    )[signature_name]

    count = 0

    for x, y in orbit:
        sig = signatures(x, y, m)[signature_name]

        if sig == actual_sig:
            count += 1

    fraction = count / len(orbit)

    return actual_sig, count, len(orbit), fraction


# ---------------------------------------------------------------------------
# Exact invariance test
# ---------------------------------------------------------------------------

def exact_invariance(orbit, m, signature_name):
    """
    Does the signature take exactly one value across the entire product orbit?
    """

    values = set()

    for x, y in orbit:
        values.add(signatures(x, y, m)[signature_name])

        if len(values) > 1:
            return False, len(values)

    return True, 1


# ---------------------------------------------------------------------------
# Orbit analysis for one target
# ---------------------------------------------------------------------------

def analyze_target_modulus(p, q, m, target_index=None, verbose=False):
    n_mod = (p * q) % m
    actual = target_residue_pair(p, q, m)

    if gcd(n_mod, m) != 1:
        return None

    orbit = product_orbit(n_mod, m)

    if actual[0] * actual[1] % m != n_mod:
        raise AssertionError("actual target is not in product orbit")

    actual_in_orbit = actual in orbit

    if not actual_in_orbit:
        # Orientation can differ, but for ordered orbit it should still be there.
        raise AssertionError("actual factor pair missing from orbit")

    sig_names = [
        "sum",
        "diff2",
        "cube_sum",
        "cube_prod",
        "Fsum",
        "Fprod",
        "Fdiff",
        "xp1_yp1",
        "xp1_sum",
        "gcd_xp1",
        "gcd_yp1",
        "order_pair",
        "order_product",
        "order_gcd",
        "order_lcm",
        "cube_class_pair",
        "F_zero_pair",
        "Jacobi_pair",
    ]

    rows = []

    for name in sig_names:
        actual_sig, count, total, fraction = signature_rarity(
            orbit,
            actual,
            m,
            name
        )

        invariant, distinct = exact_invariance(
            orbit,
            m,
            name
        )

        rows.append({
            "name": name,
            "actual": actual_sig,
            "count": count,
            "total": total,
            "fraction": fraction,
            "invariant": invariant,
            "distinct": distinct,
        })

    if verbose:
        print()
        print("TARGET", target_index)
        print("  p =", p)
        print("  q =", q)
        print("  n mod m =", n_mod)
        print("  actual pair =", actual)
        print("  orbit size =", len(orbit))

        for row in rows:
            print(
                f"  {row['name']:18s}"
                f" actual={str(row['actual']):24s}"
                f" count={row['count']:5d}/{row['total']:<5d}"
                f" fraction={row['fraction']:.6f}"
                f" distinct={row['distinct']}"
            )

    return {
        "n_mod": n_mod,
        "actual": actual,
        "orbit": orbit,
        "rows": rows,
    }


# ---------------------------------------------------------------------------
# Aggregate target statistics
# ---------------------------------------------------------------------------

def aggregate_modulus_results(targets, m):
    """
    Analyze all targets and summarize how rare the actual factor-side
    signature is within its complete product orbit.
    """

    valid = 0
    all_rows = defaultdict(list)

    for idx, (p, q) in enumerate(targets, start=1):
        result = analyze_target_modulus(
            p, q, m,
            target_index=idx,
            verbose=False
        )

        if result is None:
            continue

        valid += 1

        for row in result["rows"]:
            all_rows[row["name"]].append(row)

    return valid, all_rows


# ---------------------------------------------------------------------------
# Combined modulus
# ---------------------------------------------------------------------------

def combined_modulus(r_list):
    M = 1

    for r in r_list:
        M *= F(r)

    return M


def combined_signature(x, y, r_list):
    """
    Signature across several F(r) moduli.

    This deliberately contains both n-derived and factor-side quantities.
    """

    result = []

    for r in r_list:
        m = F(r)

        result.append((
            x % m,
            y % m,
            (x + y) % m,
            (x - y) % m,
            pow(x, 3, m),
            pow(y, 3, m),
            cyclotomic_F(x) % m,
            cyclotomic_F(y) % m,
            multiplicative_order(x % m, m),
            multiplicative_order(y % m, m),
        ))

    return tuple(result)


def combined_product_orbit(p, q, r_list):
    """
    Enumerate the complete orbit modulo the combined modulus.

    Only feasible when M is small enough.
    """

    M = combined_modulus(r_list)

    if M > MAX_EXHAUSTIVE_M:
        return None

    n_mod = (p * q) % M

    if gcd(n_mod, M) != 1:
        return None

    U = units(M)

    orbit = []

    for x in U:
        y = (n_mod * pow(x, -1, M)) % M
        orbit.append((x, y))

    return orbit


# ---------------------------------------------------------------------------
# Combined orbit analysis
# ---------------------------------------------------------------------------

def analyze_combined_stage(targets, r_list):
    M = combined_modulus(r_list)

    print()
    print("=" * 78)
    print("COMBINED ORBIT STAGE")
    print("=" * 78)
    print("r values =", r_list)
    print("F values =", [F(r) for r in r_list])
    print("M =", M)
    print("M bits =", M.bit_length())

    if M > MAX_EXHAUSTIVE_M:
        print("EXHAUSTIVE ORBIT = SKIPPED")
        print("reason: M exceeds MAX_EXHAUSTIVE_M")
        return

    U = units(M)

    print("unit count =", len(U))

    signature_rarities = []

    valid_targets = 0

    for idx, (p, q) in enumerate(targets, start=1):

        n_mod = (p * q) % M

        if gcd(n_mod, M) != 1:
            print(
                f"TARGET {idx:2d}: skipped "
                f"(n not a unit modulo M)"
            )
            continue

        valid_targets += 1

        actual = (p % M, q % M)

        orbit = combined_product_orbit(
            p, q, r_list
        )

        actual_sig = combined_signature(
            actual[0],
            actual[1],
            r_list
        )

        matching = 0

        for x, y in orbit:
            if combined_signature(x, y, r_list) == actual_sig:
                matching += 1

        fraction = matching / len(orbit)

        signature_rarities.append(fraction)

        print(
            f"TARGET {idx:2d}"
            f" orbit={len(orbit):5d}"
            f" actual_signature_count={matching:5d}"
            f" fraction={fraction:.8f}"
        )

    if signature_rarities:
        print()
        print("COMBINED-STAGE SUMMARY")
        print("valid targets =", valid_targets)
        print(
            "min actual-signature fraction =",
            f"{min(signature_rarities):.8f}"
        )
        print(
            "max actual-signature fraction =",
            f"{max(signature_rarities):.8f}"
        )
        print(
            "avg actual-signature fraction =",
            f"{sum(signature_rarities)/len(signature_rarities):.8f}"
        )


# ---------------------------------------------------------------------------
# Adversarial witness search
# ---------------------------------------------------------------------------

def find_adversarial_witness(p, q, m, signature_name):
    """
    Find an orbit pair having the same n residue but a different signature
    from the actual factor pair.
    """

    n_mod = (p * q) % m

    if gcd(n_mod, m) != 1:
        return None

    actual = (p % m, q % m)

    actual_sig = signatures(
        actual[0],
        actual[1],
        m
    )[signature_name]

    orbit = product_orbit(n_mod, m)

    for x, y in orbit:
        sig = signatures(x, y, m)[signature_name]

        if sig != actual_sig:
            return {
                "actual": actual,
                "alternative": (x, y),
                "n_mod": n_mod,
                "actual_signature": actual_sig,
                "alternative_signature": sig,
            }

    return None


# ---------------------------------------------------------------------------
# Special cyclotomic statistics
# ---------------------------------------------------------------------------

def cyclotomic_statistics(targets, m):
    """
    Compare actual factors against complete product orbits for several
    particularly interesting quantities.
    """

    names = [
        "sum",
        "cube_sum",
        "Fsum",
        "Fprod",
        "Fdiff",
        "order_pair",
        "order_gcd",
        "order_lcm",
        "F_zero_pair",
        "Jacobi_pair",
    ]

    print()
    print("=" * 78)
    print("ACTUAL FACTOR RARITY AGAINST COMPLETE PRODUCT ORBIT")
    print("=" * 78)
    print("m =", m)

    aggregate = {
        name: []
        for name in names
    }

    for idx, (p, q) in enumerate(targets, start=1):

        result = analyze_target_modulus(
            p, q, m,
            target_index=idx,
            verbose=False
        )

        if result is None:
            continue

        lookup = {
            row["name"]: row
            for row in result["rows"]
        }

        for name in names:
            aggregate[name].append(
                lookup[name]["fraction"]
            )

    print()
    print(
        f"{'signature':18s}"
        f"{'min':>12s}"
        f"{'median':>12s}"
        f"{'mean':>12s}"
        f"{'max':>12s}"
    )

    print("-" * 66)

    for name in names:
        vals = sorted(aggregate[name])

        if not vals:
            continue

        mid = len(vals) // 2

        if len(vals) % 2:
            median = vals[mid]
        else:
            median = (vals[mid - 1] + vals[mid]) / 2

        print(
            f"{name:18s}"
            f"{min(vals):12.6f}"
            f"{median:12.6f}"
            f"{sum(vals)/len(vals):12.6f}"
            f"{max(vals):12.6f}"
        )


# ---------------------------------------------------------------------------
# Exact invariant classification
# ---------------------------------------------------------------------------

def invariant_classification(m):
    """
    Test all candidate quantities against every product orbit in the
    complete unit group.

    This answers:

        Is the quantity genuinely a function of xy mod m?

    without relying on random samples.
    """

    print()
    print("=" * 78)
    print("EXACT UNIT-ACTION INVARIANT CLASSIFICATION")
    print("=" * 78)
    print("m =", m)

    U = units(m)

    if m > MAX_EXHAUSTIVE_M:
        print("SKIPPED: modulus exceeds exhaustive limit")
        return

    signature_names = [
        "sum",
        "diff2",
        "cube_sum",
        "cube_prod",
        "Fsum",
        "Fprod",
        "Fdiff",
        "xp1_yp1",
        "xp1_sum",
        "gcd_xp1",
        "gcd_yp1",
        "order_pair",
        "order_product",
        "order_gcd",
        "order_lcm",
        "cube_class_pair",
        "F_zero_pair",
        "Jacobi_pair",
    ]

    # Group all ordered unit pairs by their product.
    product_classes = defaultdict(list)

    for x in U:
        for y in U:
            product_classes[(x * y) % m].append((x, y))

    print("unit count =", len(U))
    print("product classes =", len(product_classes))

    for name in signature_names:

        invariant_classes = 0
        total_classes = 0
        maximum_values = 0

        for product, pairs in product_classes.items():

            values = set()

            for x, y in pairs:
                values.add(
                    signatures(x, y, m)[name]
                )

            total_classes += 1

            if len(values) == 1:
                invariant_classes += 1

            maximum_values = max(
                maximum_values,
                len(values)
            )

        fraction = invariant_classes / total_classes

        print(
            f"{name:18s}"
            f" invariant={invariant_classes:5d}/{total_classes:<5d}"
            f" fraction={fraction:.6f}"
            f" max_values_per_product={maximum_values}"
        )


# ---------------------------------------------------------------------------
# Target collision analysis
# ---------------------------------------------------------------------------

def target_collision_analysis(targets, r_list):
    """
    Reproduce the earlier target-group effect, but then compare it against
    complete orbit ambiguity.

    The key question:

        Does a signature uniquely identify the target because it is
        genuinely factor-restrictive, or simply because the sample is small?
    """

    M = combined_modulus(r_list)

    print()
    print("=" * 78)
    print("TARGET-LEVEL COLLISION VS ORBIT AMBIGUITY")
    print("=" * 78)
    print("r values =", r_list)
    print("M =", M)

    target_groups = defaultdict(list)

    valid = []

    for idx, (p, q) in enumerate(targets, start=1):

        n_mod = (p * q) % M

        if gcd(n_mod, M) != 1:
            continue

        pair = (p % M, q % M)

        sig = combined_signature(
            pair[0],
            pair[1],
            r_list
        )

        target_groups[sig].append(idx)
        valid.append((idx, p, q, n_mod, sig))

    print("valid targets =", len(valid))
    print("distinct signatures =", len(target_groups))

    ambiguous = [
        group
        for group in target_groups.values()
        if len(group) > 1
    ]

    print("ambiguous target signatures =", len(ambiguous))

    if ambiguous:
        print("examples:")
        for group in ambiguous[:8]:
            print(" ", group)

    # Now show orbit-level uniqueness.
    print()
    print("ORBIT-LEVEL TEST")

    for idx, p, q, n_mod, actual_sig in valid:

        orbit = combined_product_orbit(
            p, q, r_list
        )

        if orbit is None:
            continue

        count = 0

        for x, y in orbit:
            if combined_signature(x, y, r_list) == actual_sig:
                count += 1

        print(
            f"target={idx:2d}"
            f" orbit={len(orbit):5d}"
            f" same-signature={count:5d}"
            f" fraction={count/len(orbit):.8f}"
        )


# ---------------------------------------------------------------------------
# Print representative witness
# ---------------------------------------------------------------------------

def print_representative_witness(targets):
    p, q = targets[0]

    print()
    print("=" * 78)
    print("REPRESENTATIVE TARGET")
    print("=" * 78)

    n = p * q
    s = p + q
    delta = (p - q) ** 2

    print("p     =", p)
    print("q     =", q)
    print("n     =", n)
    print("s     =", s)
    print("Delta =", delta)
    print("n bits =", n.bit_length())


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("KAPPA NEXT EXPERIMENT")
    print("FULL ORBIT VS ACTUAL FACTOR DISTRIBUTION")
    print("ADVERSARIAL CYCLOTOMIC SIGNATURE TEST")
    print("=" * 78)

    print("random seed =", SEED)
    print("R values    =", R_VALUES)
    print("targets     =", TARGET_COUNT)
    print("CSV output  = NONE")
    print()

    targets = generate_targets(TARGET_COUNT)

    print_representative_witness(targets)

    # -----------------------------------------------------------------------
    # 1. Modulus inventory
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. MODULUS INVENTORY")
    print("=" * 78)

    for r in R_VALUES:
        m = F(r)

        U = units(m)

        roots_F = [
            x for x in range(m)
            if cyclotomic_F(x) % m == 0
        ]

        roots_cube_minus_one = [
            x for x in range(m)
            if pow(x, 3, m) == (m - 1) % m
        ]

        print(
            f"r={r:2d}"
            f" m={m:6d}"
            f" units={len(U):6d}"
            f" F-roots={len(roots_F):4d}"
            f" x^3=-1 roots={len(roots_cube_minus_one):4d}"
        )

    # -----------------------------------------------------------------------
    # 2. Exact invariant classification
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT PRODUCT-ORBIT INVARIANT SEARCH")
    print("=" * 78)

    # Do the smaller moduli exhaustively.
    for r in R_VALUES:
        m = F(r)

        if m <= MAX_EXHAUSTIVE_M:
            invariant_classification(m)

    # -----------------------------------------------------------------------
    # 3. Actual factor rarity
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. ACTUAL FACTOR VS COMPLETE PRODUCT ORBIT")
    print("=" * 78)

    selected_r = [
        2, 3, 5, 7, 11, 13, 17, 19, 23
    ]

    for r in selected_r:
        m = F(r)
        cyclotomic_statistics(targets, m)

    # -----------------------------------------------------------------------
    # 4. Explicit adversarial witnesses
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. ADVERSARIAL WITNESS SEARCH")
    print("=" * 78)

    witness_names = [
        "sum",
        "cube_sum",
        "Fsum",
        "Fprod",
        "Fdiff",
        "order_pair",
        "order_gcd",
        "order_lcm",
    ]

    for r in [2, 3, 5, 7, 11]:

        m = F(r)

        print()
        print(f"r={r} m={m}")

        for name in witness_names:

            witness = find_adversarial_witness(
                targets[0][0],
                targets[0][1],
                m,
                name
            )

            if witness is None:
                print(
                    f"  {name:14s}: invariant "
                    f"(no adversarial witness)"
                )
            else:
                print(
                    f"  {name:14s}: "
                    f"actual={witness['actual']}"
                    f" alt={witness['alternative']}"
                    f" actual_sig={witness['actual_signature']}"
                    f" alt_sig={witness['alternative_signature']}"
                )

    # -----------------------------------------------------------------------
    # 5. Combined modulus stages
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. COMBINED-MODULUS ORBIT TEST")
    print("=" * 78)

    for stage in COMBINED_STAGES:

        M = combined_modulus(stage)

        # Only exhaustive stages.
        if M <= MAX_EXHAUSTIVE_M:
            analyze_combined_stage(
                targets,
                stage
            )
        else:
            print()
            print(
                "COMBINED STAGE",
                stage,
                "M=",
                M,
                "SKIPPED: too large for exhaustive orbit"
            )

    # -----------------------------------------------------------------------
    # 6. Re-test the apparent target collapse
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. TARGET COLLISION VS ORBIT COLLISION")
    print("=" * 78)

    for stage in COMBINED_STAGES:

        M = combined_modulus(stage)

        if M <= MAX_EXHAUSTIVE_M:
            target_collision_analysis(
                targets,
                stage
            )

    # -----------------------------------------------------------------------
    # 7. Special F(p)=0 statistics
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. SPECIAL EVENT: F(p)=0 OR F(q)=0")
    print("=" * 78)

    for r in R_VALUES:

        m = F(r)

        p_hits = 0
        q_hits = 0
        both_hits = 0

        for p, q in targets:

            ph = F_zero_flag(p % m, m)
            qh = F_zero_flag(q % m, m)

            if ph:
                p_hits += 1

            if qh:
                q_hits += 1

            if ph and qh:
                both_hits += 1

        roots = sum(
            1 for x in range(m)
            if cyclotomic_F(x) % m == 0
        )

        print(
            f"r={r:2d}"
            f" m={m:6d}"
            f" roots={roots:3d}"
            f" p_hits={p_hits:2d}/{len(targets)}"
            f" q_hits={q_hits:2d}/{len(targets)}"
            f" both={both_hits:2d}"
        )

    # -----------------------------------------------------------------------
    # 8. Final interpretation
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. FINAL CLASSIFICATION")
    print("=" * 78)

    print()
    print("The decisive statistic in this experiment is:")
    print()
    print("    actual-signature orbit fraction")
    print()
    print("For an actual target n=p*q, this is")
    print()
    print("    # admissible factor pairs sharing the actual signature")
    print("    ------------------------------------------------------")
    print("                    # all factor pairs with xy=n")
    print()
    print("Interpretation:")
    print()
    print("  FRACTION ~ 1")
    print("      Signature is effectively n-derived / non-restrictive.")
    print()
    print("  FRACTION small")
    print("      Actual factors occupy a restricted subset of the product orbit.")
    print("      This is the interesting case.")
    print()
    print("  EXACTLY ONE")
    print("      The signature uniquely identifies the factor residue pair")
    print("      inside the complete product orbit.")
    print()
    print("IMPORTANT:")
    print("A signature being unique among the 24 sampled targets is NOT enough.")
    print("It must also be unusually selective inside the complete adversarial")
    print("product orbit.")
    print()
    print("If all nontrivial signatures have large orbit fractions, then")
    print("cyclotomic residue classes are behaving like generic auxiliary")
    print("moduli and the next direction should move toward norms, resultants,")
    print("or higher cyclotomic structure.")
    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()

