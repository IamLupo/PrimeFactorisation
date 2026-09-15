#!/usr/bin/env python3

import math
import random
from collections import Counter, defaultdict

# ==============================================================================
# KAPPA NEXT EXPERIMENT
# ACTUAL FACTOR ORBIT RARITY
# COMPLETE UNIT ORBIT VS ACTUAL FACTOR PAIR
# ==============================================================================

SEED = 20260814
random.seed(SEED)

R_VALUES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

# Keep this consistent with the previous experiments.
TARGETS = 24

# Number of signatures retained in combined-modulus tests.
MAX_COMBINED_SIGNATURES = 8

# Maximum number of targets printed in detail.
DETAIL_TARGETS = 6

# ==============================================================================
# TARGET GENERATION
# ==============================================================================

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


def make_targets(count):
    """
    Deterministic semiprime targets.

    The first target is the representative target from the previous runs.
    The remaining targets are generated deterministically from the seed.
    """

    targets = [
        (3622631, 3922519)
    ]

    rng = random.Random(SEED)

    while len(targets) < count:
        p = rng.randrange(2_000_000, 5_000_000)
        q = rng.randrange(2_000_000, 5_000_000)

        if p == q:
            continue

        if not is_prime(p) or not is_prime(q):
            continue

        if p > q:
            p, q = q, p

        if (p, q) not in targets:
            targets.append((p, q))

    return targets


# ==============================================================================
# BASIC NUMBER THEORY
# ==============================================================================

def units_mod(m):
    return [x for x in range(1, m) if math.gcd(x, m) == 1]


def inv_mod(x, m):
    return pow(x, -1, m)


def multiplicative_order(x, m):
    if math.gcd(x, m) != 1:
        return None

    value = 1
    k = 0

    while True:
        k += 1
        value = (value * x) % m
        if value == 1:
            return k


def jacobi_symbol(a, n):
    """
    Jacobi symbol (a/n), n odd positive.
    """

    if n <= 0 or n % 2 == 0:
        return None

    a %= n
    result = 1

    while a:
        while a % 2 == 0:
            a //= 2
            r = n % 8
            if r in (3, 5):
                result = -result

        a, n = n, a

        if a % 4 == 3 and n % 4 == 3:
            result = -result

        a %= n

    return result if n == 1 else 0


# ==============================================================================
# CYCLOTOMIC-STYLE FUNCTION
# ==============================================================================

def cyclotomic_F(x, r, m):
    """
    F_r(x) = 1 + x + ... + x^(r-1) mod m

    Since r is prime in R_VALUES, this is the usual Phi_r(x).
    """

    x %= m

    total = 0
    power = 1

    for _ in range(r):
        total = (total + power) % m
        power = (power * x) % m

    return total


# ==============================================================================
# CUBIC CLASS
# ==============================================================================

def cubic_class(x, m):
    """
    Classifies x according to x^((phi(m))/gcd(3,phi(m))).

    This is deliberately only used for units.
    """

    if math.gcd(x, m) != 1:
        return None

    phi = len(units_mod(m))

    e = phi // math.gcd(3, phi)

    return pow(x, e, m)


# ==============================================================================
# MODULUS DEFINITION
# ==============================================================================

def modulus_for_r(r):
    """
    Previous experiment uses m = 2*r^2 - 1.
    """

    return 2 * r * r - 1


# ==============================================================================
# SIGNATURES
# ==============================================================================

def signature_values(x, y, r, m):
    Fx = cyclotomic_F(x, r, m)
    Fy = cyclotomic_F(y, r, m)

    ox = multiplicative_order(x, m)
    oy = multiplicative_order(y, m)

    cx = cubic_class(x, m)
    cy = cubic_class(y, m)

    # All signatures are deliberately factor-side.
    #
    # They are NOT assumed to be n-only invariants.
    #
    # The purpose here is to measure whether the actual factor point
    # occupies an unusually rare location in its complete product orbit.

    sig = {}

    sig["Fpair"] = (Fx, Fy)
    sig["Fsorted"] = tuple(sorted((Fx, Fy)))

    sig["Fzero_pair"] = (int(Fx == 0), int(Fy == 0))

    sig["order_pair"] = (ox, oy)
    sig["order_sorted"] = tuple(sorted((ox, oy)))

    if ox is not None and oy is not None:
        sig["order_gcd"] = math.gcd(ox, oy)
        sig["order_lcm"] = math.lcm(ox, oy)
        sig["order_product"] = ox * oy
    else:
        sig["order_gcd"] = None
        sig["order_lcm"] = None
        sig["order_product"] = None

    sig["cube_class_pair"] = (cx, cy)
    sig["cube_class_sorted"] = tuple(sorted((cx, cy)))

    sig["Jacobi_pair"] = (
        jacobi_symbol(x, m) if m % 2 else None,
        jacobi_symbol(y, m) if m % 2 else None,
    )

    sig["Fsum"] = (Fx + Fy) % m
    sig["Fprod"] = (Fx * Fy) % m
    sig["Fdiff"] = (Fx - Fy) % m

    sig["cube_sum"] = (pow(x, 3, m) + pow(y, 3, m)) % m
    sig["cube_prod"] = (pow(x, 3, m) * pow(y, 3, m)) % m

    sig["xp1_yp1"] = ((x + 1) * (y + 1)) % m
    sig["xp1_sum"] = (x + 1 + y + 1) % m

    sig["gcd_xp1"] = math.gcd(x + 1, m)
    sig["gcd_yp1"] = math.gcd(y + 1, m)

    return sig


SIGNATURE_NAMES = [
    "Fpair",
    "Fsorted",
    "Fzero_pair",
    "order_pair",
    "order_sorted",
    "order_gcd",
    "order_lcm",
    "order_product",
    "cube_class_pair",
    "cube_class_sorted",
    "Jacobi_pair",
    "Fsum",
    "Fprod",
    "Fdiff",
    "cube_sum",
    "cube_prod",
    "xp1_yp1",
    "xp1_sum",
    "gcd_xp1",
    "gcd_yp1",
]


# ==============================================================================
# COMPLETE PRODUCT ORBIT
# ==============================================================================

def build_orbit(p, q, m, units):
    """
    Complete unit orbit:

        (p*t, q*t^-1) mod m

    This preserves p*q mod m exactly.
    """

    px = p % m
    qx = q % m

    orbit = []

    for t in units:
        x = (px * t) % m
        y = (qx * inv_mod(t, m)) % m

        orbit.append((x, y))

    return orbit


# ==============================================================================
# ORBIT SIGNATURE DISTRIBUTION
# ==============================================================================

def orbit_distributions(orbit, r, m):
    distributions = {
        name: Counter()
        for name in SIGNATURE_NAMES
    }

    signature_cache = []

    for x, y in orbit:
        sig = signature_values(x, y, r, m)
        signature_cache.append(sig)

        for name in SIGNATURE_NAMES:
            distributions[name][sig[name]] += 1

    return distributions, signature_cache


# ==============================================================================
# ACTUAL RARITY
# ==============================================================================

def rarity_for_signature(counter, actual_value, orbit_size):
    hits = counter.get(actual_value, 0)

    if orbit_size == 0:
        return 0, 0.0

    return hits, hits / orbit_size


def empirical_percentile(counter, actual_value, orbit_size):
    """
    Lower percentile means the actual signature value is in a small orbit
    class.

    For non-ordered signatures, we measure only class frequency.

    The returned quantity is:

        fraction of orbit points belonging to the actual signature class.
    """

    hits = counter.get(actual_value, 0)

    if orbit_size == 0:
        return 0.0

    return hits / orbit_size


# ==============================================================================
# ONE MODULUS ANALYSIS
# ==============================================================================

def analyze_modulus(p, q, r, verbose=False):
    m = modulus_for_r(r)

    units = units_mod(m)

    # Actual residues.
    x_actual = p % m
    y_actual = q % m

    # Complete orbit.
    orbit = build_orbit(p, q, m, units)

    # Sanity check: product must remain constant.
    target_product = (x_actual * y_actual) % m

    for x, y in orbit:
        if (x * y) % m != target_product:
            raise RuntimeError(
                f"Orbit product failure: r={r}, m={m}, "
                f"x={x}, y={y}"
            )

    distributions, signature_cache = orbit_distributions(
        orbit, r, m
    )

    actual_sig = signature_values(
        x_actual,
        y_actual,
        r,
        m
    )

    result = {
        "r": r,
        "m": m,
        "orbit_size": len(orbit),
        "units": len(units),
        "x": x_actual,
        "y": y_actual,
        "actual": {},
    }

    for name in SIGNATURE_NAMES:
        hits, fraction = rarity_for_signature(
            distributions[name],
            actual_sig[name],
            len(orbit)
        )

        result["actual"][name] = {
            "value": actual_sig[name],
            "hits": hits,
            "fraction": fraction,
        }

    if verbose:
        print()
        print("-" * 78)
        print(f"r={r:2d}  m={m:5d}  orbit={len(orbit):5d}")
        print(
            f"actual residues: "
            f"x={x_actual:5d} y={y_actual:5d} "
            f"xy={target_product:5d}"
        )
        print("-" * 78)

        for name in SIGNATURE_NAMES:
            row = result["actual"][name]

            print(
                f"{name:20s} "
                f"hits={row['hits']:5d}/{len(orbit):5d} "
                f"freq={row['fraction']:.6f} "
                f"value={row['value']}"
            )

    return result


# ==============================================================================
# RANK / RARITY SUMMARY
# ==============================================================================

def print_rarity_table(result, names):
    print()
    print(
        f"r={result['r']:2d} m={result['m']:5d} "
        f"orbit={result['orbit_size']:5d}"
    )

    print(
        f"{'signature':22s}"
        f"{'hits':>10s}"
        f"{'frequency':>13s}"
        f"{'actual value':>24s}"
    )

    print("-" * 72)

    for name in names:
        row = result["actual"][name]

        print(
            f"{name:22s}"
            f"{row['hits']:>5d}/{result['orbit_size']:<4d}"
            f"{row['fraction']:>13.6f}"
            f"{str(row['value']):>24s}"
        )


# ==============================================================================
# CROSS-MODULUS COMBINATION
# ==============================================================================

def combined_actual_rarity(
    p,
    q,
    selected_r,
    signature_name
):
    """
    Combine signatures across several moduli.

    For each unit-orbit point at each modulus, record the actual signature
    class frequency.

    We deliberately use the product of marginal frequencies as a diagnostic,
    NOT as an assertion of independence.
    """

    values = []

    for r in selected_r:
        result = analyze_modulus(
            p, q, r, verbose=False
        )

        row = result["actual"][signature_name]

        values.append({
            "r": r,
            "m": result["m"],
            "hits": row["hits"],
            "orbit": result["orbit_size"],
            "frequency": row["fraction"],
        })

    product_frequency = 1.0

    for item in values:
        product_frequency *= item["frequency"]

    return values, product_frequency


# ==============================================================================
# TARGET-LEVEL ANALYSIS
# ==============================================================================

def analyze_target(target_index, p, q):
    n = p * q
    s = p + q
    delta = (p - q) ** 2

    print()
    print("=" * 78)
    print(f"TARGET {target_index}")
    print("=" * 78)
    print(f"p       = {p}")
    print(f"q       = {q}")
    print(f"n       = {n}")
    print(f"s       = {s}")
    print(f"Delta   = {delta}")
    print(f"n bits  = {n.bit_length()}")

    # ------------------------------------------------------------
    # Per-modulus analysis
    # ------------------------------------------------------------

    results = []

    for r in R_VALUES:
        result = analyze_modulus(
            p,
            q,
            r,
            verbose=False
        )

        results.append(result)

    # ------------------------------------------------------------
    # Print selected detailed moduli.
    # ------------------------------------------------------------

    print()
    print("PER-MODULUS ACTUAL ORBIT RARITY")
    print("-" * 78)

    selected_names = [
        "Fpair",
        "Fsorted",
        "Fzero_pair",
        "order_pair",
        "order_sorted",
        "order_lcm",
        "cube_class_pair",
        "Jacobi_pair",
        "Fsum",
        "Fprod",
        "Fdiff",
        "cube_sum",
        "cube_prod",
    ]

    for result in results:
        print_rarity_table(
            result,
            selected_names
        )

    # ------------------------------------------------------------
    # Rarity counts.
    # ------------------------------------------------------------

    print()
    print("=" * 78)
    print("RARITY COUNTS ACROSS MODULI")
    print("=" * 78)

    print(
        f"{'signature':22s}"
        f"{'<1%':>8s}"
        f"{'<5%':>8s}"
        f"{'<10%':>9s}"
        f"{'median':>12s}"
    )

    print("-" * 65)

    for name in selected_names:
        fractions = [
            r["actual"][name]["fraction"]
            for r in results
        ]

        less_1 = sum(x < 0.01 for x in fractions)
        less_5 = sum(x < 0.05 for x in fractions)
        less_10 = sum(x < 0.10 for x in fractions)

        ordered = sorted(fractions)
        mid = len(ordered) // 2

        if len(ordered) % 2:
            median = ordered[mid]
        else:
            median = (
                ordered[mid - 1] + ordered[mid]
            ) / 2

        print(
            f"{name:22s}"
            f"{less_1:>8d}"
            f"{less_5:>8d}"
            f"{less_10:>9d}"
            f"{median:>12.6f}"
        )

    # ------------------------------------------------------------
    # Combined modulus tests.
    # ------------------------------------------------------------

    print()
    print("=" * 78)
    print("CROSS-MODULUS ACTUAL RARITY")
    print("=" * 78)

    combinations = [
        [2, 3],
        [2, 3, 5],
        [2, 3, 5, 7],
        [2, 3, 5, 7, 11],
    ]

    combined_names = [
        "Fpair",
        "Fsorted",
        "order_pair",
        "order_sorted",
        "order_lcm",
        "cube_class_pair",
        "Jacobi_pair",
        "Fsum",
    ]

    print(
        f"{'signature':22s}"
        f"{'2,3':>14s}"
        f"{'2,3,5':>14s}"
        f"{'2,3,5,7':>16s}"
        f"{'2,3,5,7,11':>20s}"
    )

    print("-" * 90)

    for name in combined_names:
        outputs = []

        for rs in combinations:
            vals, product_frequency = combined_actual_rarity(
                p,
                q,
                rs,
                name
            )

            outputs.append(product_frequency)

        print(
            f"{name:22s}"
            f"{outputs[0]:14.8g}"
            f"{outputs[1]:14.8g}"
            f"{outputs[2]:16.8g}"
            f"{outputs[3]:20.8g}"
        )

    return results


# ==============================================================================
# GLOBAL SUMMARY
# ==============================================================================

def global_summary(all_target_results):
    print()
    print("=" * 78)
    print("GLOBAL ACTUAL-FACTOR RARITY SUMMARY")
    print("=" * 78)

    names = [
        "Fpair",
        "Fsorted",
        "Fzero_pair",
        "order_pair",
        "order_sorted",
        "order_lcm",
        "cube_class_pair",
        "Jacobi_pair",
        "Fsum",
        "Fprod",
        "Fdiff",
        "cube_sum",
        "cube_prod",
    ]

    print(
        f"{'signature':22s}"
        f"{'target-moduli':>14s}"
        f"{'<1%':>8s}"
        f"{'<5%':>8s}"
        f"{'<10%':>9s}"
        f"{'median freq':>14s}"
    )

    print("-" * 82)

    for name in names:
        frequencies = []

        for target_results in all_target_results:
            for result in target_results:
                frequencies.append(
                    result["actual"][name]["fraction"]
                )

        if not frequencies:
            continue

        less_1 = sum(x < 0.01 for x in frequencies)
        less_5 = sum(x < 0.05 for x in frequencies)
        less_10 = sum(x < 0.10 for x in frequencies)

        ordered = sorted(frequencies)
        mid = len(ordered) // 2

        if len(ordered) % 2:
            median = ordered[mid]
        else:
            median = (
                ordered[mid - 1] +
                ordered[mid]
            ) / 2

        print(
            f"{name:22s}"
            f"{len(frequencies):>14d}"
            f"{less_1:>8d}"
            f"{less_5:>8d}"
            f"{less_10:>9d}"
            f"{median:>14.6f}"
        )


# ==============================================================================
# SPECIAL TEST:
# ACTUAL FACTOR VS RANDOM PRODUCT-COMPATIBLE POINT
# ==============================================================================

def random_orbit_baseline(
    p,
    q,
    r,
    trials=500
):
    """
    Draw random unit-action points from the same product orbit.

    This gives a direct Monte-Carlo baseline.

    The actual point is compared against randomly selected points from
    precisely the same orbit, so the baseline does not change n mod m.
    """

    m = modulus_for_r(r)
    units = units_mod(m)

    px = p % m
    qx = q % m

    actual_x = px
    actual_y = qx

    actual_sig = signature_values(
        actual_x,
        actual_y,
        r,
        m
    )

    rng = random.Random(
        SEED + r + (p % 1000003)
    )

    random_frequencies = defaultdict(list)

    for _ in range(trials):
        t = rng.choice(units)

        x = (px * t) % m
        y = (qx * inv_mod(t, m)) % m

        sig = signature_values(x, y, r, m)

        for name in SIGNATURE_NAMES:
            random_frequencies[name].append(
                sig[name]
            )

    print()
    print("-" * 78)
    print(
        f"RANDOM SAME-ORBIT BASELINE "
        f"r={r} m={m} trials={trials}"
    )
    print("-" * 78)

    names = [
        "Fpair",
        "Fsorted",
        "order_pair",
        "order_sorted",
        "order_lcm",
        "cube_class_pair",
        "Jacobi_pair",
        "Fsum",
        "Fprod",
        "Fdiff",
        "cube_sum",
        "cube_prod",
    ]

    for name in names:
        samples = random_frequencies[name]

        actual = actual_sig[name]
        hits = sum(v == actual for v in samples)

        print(
            f"{name:22s}"
            f" actual={str(actual):>18s}"
            f" random_hits={hits:4d}/{trials}"
            f" freq={hits/trials:.4f}"
        )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print("KAPPA NEXT EXPERIMENT")
    print("ACTUAL FACTOR ORBIT RARITY")
    print("COMPLETE UNIT ORBIT VS ACTUAL FACTOR PAIR")
    print("=" * 78)

    print(f"random seed = {SEED}")
    print(f"R values    = {R_VALUES}")
    print(f"targets     = {TARGETS}")
    print("CSV output  = NONE")

    # ------------------------------------------------------------
    # Modulus inventory
    # ------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. MODULUS INVENTORY")
    print("=" * 78)

    for r in R_VALUES:
        m = modulus_for_r(r)
        units = units_mod(m)

        print(
            f"r={r:2d} "
            f"m={m:6d} "
            f"units={len(units):6d}"
        )

    # ------------------------------------------------------------
    # Targets
    # ------------------------------------------------------------

    targets = make_targets(TARGETS)

    all_target_results = []

    # ------------------------------------------------------------
    # Detailed target analysis
    # ------------------------------------------------------------

    for i, (p, q) in enumerate(targets, start=1):

        # Full analysis for all targets.
        #
        # Only the first DETAIL_TARGETS are printed in complete detail.
        #
        # The calculations are still performed for all targets.
        if i <= DETAIL_TARGETS:
            results = analyze_target(i, p, q)
        else:
            results = []

            for r in R_VALUES:
                results.append(
                    analyze_modulus(
                        p,
                        q,
                        r,
                        verbose=False
                    )
                )

        all_target_results.append(results)

        # --------------------------------------------------------
        # Lightweight target summary for later targets.
        # --------------------------------------------------------

        if i > DETAIL_TARGETS:

            print()
            print(
                f"TARGET {i:2d} "
                f"n={p*q} "
                f"bits={(p*q).bit_length()}"
            )

            for name in [
                "Fpair",
                "Fsorted",
                "order_pair",
                "order_lcm",
                "cube_class_pair",
            ]:

                frequencies = [
                    x["actual"][name]["fraction"]
                    for x in results
                ]

                less_1 = sum(
                    x < 0.01
                    for x in frequencies
                )

                less_5 = sum(
                    x < 0.05
                    for x in frequencies
                )

                ordered = sorted(frequencies)

                mid = len(ordered) // 2

                if len(ordered) % 2:
                    median = ordered[mid]
                else:
                    median = (
                        ordered[mid - 1]
                        + ordered[mid]
                    ) / 2

                print(
                    f"  {name:20s}"
                    f" <1%={less_1:2d}"
                    f" <5%={less_5:2d}"
                    f" median={median:.6f}"
                )

    # ------------------------------------------------------------
    # Global summary
    # ------------------------------------------------------------

    global_summary(all_target_results)

    # ------------------------------------------------------------
    # Same-orbit random baseline on representative target.
    # ------------------------------------------------------------

    p, q = targets[0]

    print()
    print("=" * 78)
    print("SAME-ORBIT RANDOM BASELINE")
    print("=" * 78)

    for r in [2, 3, 5, 7, 11, 13, 17]:
        random_orbit_baseline(
            p,
            q,
            r,
            trials=500
        )

    # ------------------------------------------------------------
    # Final interpretation
    # ------------------------------------------------------------

    print()
    print("=" * 78)
    print("INTERPRETATION GUIDE")
    print("=" * 78)

    print(
        """
The central statistic is:

    frequency =
        number of orbit points having the actual signature
        ---------------------------------------------------
                     complete orbit size

This is NOT an n-only invariant test.

It asks whether the actual factor pair is unusually rare
inside the complete unit orbit compatible with n mod m.

Interpretation:

  frequency ~ 1
      signature is common; weak evidence.

  frequency ~ 0.1
      moderately selective.

  frequency < 0.05
      potentially interesting.

  frequency < 0.01
      strongly selective at that modulus.

  repeated <1% behavior across many independent moduli
      potentially significant.

However:

  rarity alone does NOT imply factorization.

The decisive next criterion is whether the rarity persists
across independent targets and survives combined-modulus
testing without being explained by the size of the orbit.

The perfect cube_prod result is expected because it is a
direct product-derived quantity and therefore is not evidence
of factor-side information.
"""
    )

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()

