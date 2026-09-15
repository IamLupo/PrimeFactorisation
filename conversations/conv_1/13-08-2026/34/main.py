import random
import time
from collections import defaultdict

import sympy


# ==============================================================================
# KAPPA EXPERIMENT 34
# RESIDUE-CLASS SIGNATURE DIAGNOSTIC
# NO CSV OUTPUT
# ==============================================================================

SEED = 20260814

TARGET_COUNT = 12
POOL_SIZE = 5000

PRIME_MIN = 2_000_000
PRIME_MAX = 4_200_000

PREFIXES = [3, 4, 5, 6, 7]

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17,
    19, 23, 29, 31, 37, 41, 43, 47
]

random.seed(SEED)


# ==============================================================================
# PRINT
# ==============================================================================

def header(text):
    print()
    print("=" * 78)
    print(text)
    print("=" * 78)


def section(text):
    print()
    print("-" * 78)
    print(text)
    print("-" * 78)


# ==============================================================================
# MODULUS
# ==============================================================================

def build_modulus(r):
    return 3 * r + 1


# ==============================================================================
# PRIME POOL
# ==============================================================================

def generate_prime_pool():

    primes = set()

    while len(primes) < POOL_SIZE:

        x = random.randint(
            PRIME_MIN,
            PRIME_MAX
        )

        if sympy.isprime(x):
            primes.add(x)

    return sorted(primes)


# ==============================================================================
# TARGETS
# ==============================================================================

def generate_targets(prime_pool):

    targets = []
    used = set()

    while len(targets) < TARGET_COUNT:

        p, q = random.sample(
            prime_pool,
            2
        )

        if p > q:
            p, q = q, p

        pair = (p, q)

        if pair in used:
            continue

        used.add(pair)

        targets.append(
            (p, q, p * q)
        )

    return targets


# ==============================================================================
# PRIME RESIDUE GROUPS
# ==============================================================================

def build_residue_groups(prime_pool, m):

    groups = defaultdict(list)

    for p in prime_pool:
        groups[p % m].append(p)

    return groups


# ==============================================================================
# PAIR COUNT FOR TWO RESIDUE CLASSES
# ==============================================================================

def pair_count_same_class(count):
    return count * (count - 1) // 2


# ==============================================================================
# SIGNATURES ON RESIDUES
# ==============================================================================

def sig_fpair(a, b, m):
    return tuple(sorted((a, b)))


def sig_fsorted(a, b, m):
    return tuple(sorted((a, b)))


def sig_fdiff(a, b, m):
    return (a - b) % m


def sig_fsum(a, b, m):
    return (a + b) % m


def sig_fprod(a, b, m):
    return (a * b) % m


def sig_cube_sum(a, b, m):
    return (
        pow(a, 3, m) +
        pow(b, 3, m)
    ) % m


# ==============================================================================
# ORDER TABLE
# ==============================================================================

def build_order_table(m):

    table = {}

    for a in range(m):

        if sympy.gcd(a, m) != 1:
            table[a] = None
        else:
            table[a] = int(
                sympy.n_order(a, m)
            )

    return table


def sig_order_pair(a, b, m, order_table):

    oa = order_table[a]
    ob = order_table[b]

    if oa is None or ob is None:
        return None

    return tuple(sorted((oa, ob)))


# ==============================================================================
# SIGNATURE LIST
# ==============================================================================

SIGNATURES = [
    ("Fpair", sig_fpair),
    ("Fsorted", sig_fsorted),
    ("Fdiff", sig_fdiff),
    ("Fsum", sig_fsum),
    ("Fprod", sig_fprod),
    ("cube_sum", sig_cube_sum),
]


# ==============================================================================
# BUILD RESIDUE SIGNATURE COUNTS
# ==============================================================================

def build_signature_counts(
    groups,
    m,
    signature_func
):

    counts = defaultdict(int)

    residues = sorted(groups)

    for i, a in enumerate(residues):

        ca = len(groups[a])

        for j in range(i, len(residues)):

            b = residues[j]

            cb = len(groups[b])

            # Number of unordered prime pairs represented by
            # residue classes a and b.
            if a == b:
                pair_count = pair_count_same_class(ca)
            else:
                pair_count = ca * cb

            if pair_count == 0:
                continue

            sig = signature_func(
                a,
                b,
                m
            )

            counts[sig] += pair_count

    return counts


# ==============================================================================
# ORDER SIGNATURE COUNTS
# ==============================================================================

def build_order_signature_counts(
    groups,
    m,
    order_table
):

    counts = defaultdict(int)

    residues = sorted(groups)

    for i, a in enumerate(residues):

        oa = order_table[a]

        if oa is None:
            continue

        ca = len(groups[a])

        for j in range(i, len(residues)):

            b = residues[j]

            ob = order_table[b]

            if ob is None:
                continue

            cb = len(groups[b])

            if a == b:
                pair_count = pair_count_same_class(ca)
            else:
                pair_count = ca * cb

            sig = tuple(
                sorted((oa, ob))
            )

            counts[sig] += pair_count

    return counts


# ==============================================================================
# TARGET SIGNATURE
# ==============================================================================

def target_sig(
    p,
    q,
    m,
    name,
    order_table=None
):

    a = p % m
    b = q % m

    if name == "Fpair":
        return sig_fpair(a, b, m)

    if name == "Fsorted":
        return sig_fsorted(a, b, m)

    if name == "Fdiff":
        return sig_fdiff(a, b, m)

    if name == "Fsum":
        return sig_fsum(a, b, m)

    if name == "Fprod":
        return sig_fprod(a, b, m)

    if name == "cube_sum":
        return sig_cube_sum(a, b, m)

    if name == "order_pair":
        return sig_order_pair(
            a,
            b,
            m,
            order_table
        )

    raise ValueError(name)


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    start_total = time.perf_counter()

    header(
        "KAPPA EXPERIMENT 34\n"
        "RESIDUE-CLASS SIGNATURE DIAGNOSTIC\n"
        "NO CSV OUTPUT"
    )

    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGET_COUNT}")
    print(f"prime pool       = {POOL_SIZE:,}")
    print(
        f"prime interval   = "
        f"[{PRIME_MIN:,}, {PRIME_MAX:,}]"
    )
    print(f"prefixes         = {PREFIXES}")
    print(f"R values         = {R_VALUES}")

    # ==========================================================================
    # MODULI
    # ==========================================================================

    section("MODULUS INVENTORY")

    moduli = {}

    for r in R_VALUES:

        m = build_modulus(r)

        moduli[r] = m

        units = sum(
            1
            for x in range(1, m)
            if sympy.gcd(x, m) == 1
        )

        print(
            f"r={r:>3} "
            f"m={m:>7} "
            f"units={units:>7}"
        )

    # ==========================================================================
    # PRIME POOL
    # ==========================================================================

    section("PRIME POOL")

    t0 = time.perf_counter()

    prime_pool = generate_prime_pool()

    print(
        f"generated primes = {len(prime_pool):,}"
    )

    print(
        f"generation time = "
        f"{time.perf_counter() - t0:.3f}s"
    )

    # ==========================================================================
    # TARGETS
    # ==========================================================================

    section("TARGETS")

    targets = generate_targets(
        prime_pool
    )

    for i, (p, q, n) in enumerate(
        targets,
        1
    ):

        print(
            f"target {i:>2}: "
            f"p={p} "
            f"q={q} "
            f"n={n}"
        )

    # ==========================================================================
    # PRECOMPUTE RESIDUE GROUPS
    # ==========================================================================

    section(
        "RESIDUE-CLASS CACHE"
    )

    residue_groups = {}

    t0 = time.perf_counter()

    for r in R_VALUES:

        m = moduli[r]

        residue_groups[m] = (
            build_residue_groups(
                prime_pool,
                m
            )
        )

    print(
        f"residue groups ready "
        f"[{time.perf_counter() - t0:.3f}s]"
    )

    # ==========================================================================
    # ORDER CACHE
    # ==========================================================================

    section(
        "MULTIPLICATIVE ORDER TABLES"
    )

    order_tables = {}

    t0 = time.perf_counter()

    for r in R_VALUES:

        m = moduli[r]

        order_tables[m] = (
            build_order_table(m)
        )

    print(
        f"order tables ready "
        f"[{time.perf_counter() - t0:.3f}s]"
    )

    # ==========================================================================
    # CACHE SIGNATURE COUNTS
    # ==========================================================================

    section(
        "BUILDING RESIDUE SIGNATURE COUNTS"
    )

    signature_cache = {}

    for prefix in PREFIXES:

        # Prefix means first N R-values.
        selected_r = R_VALUES[:prefix]

        # For the diagnostic we test each modulus separately.
        for r in selected_r:

            m = moduli[r]

            groups = residue_groups[m]

            print(
                f"building signatures "
                f"prefix={prefix} "
                f"r={r} "
                f"m={m}"
            )

            local = {}

            t0 = time.perf_counter()

            for name, func in SIGNATURES:

                local[name] = (
                    build_signature_counts(
                        groups,
                        m,
                        func
                    )
                )

            local["order_pair"] = (
                build_order_signature_counts(
                    groups,
                    m,
                    order_tables[m]
                )
            )

            signature_cache[
                (prefix, r)
            ] = local

            print(
                f"  complete "
                f"[{time.perf_counter() - t0:.3f}s]"
            )

    # ==========================================================================
    # TARGET DIAGNOSTIC
    # ==========================================================================

    header(
        "TARGET SIGNATURE DIAGNOSTIC"
    )

    # Summary:
    #
    # wrong = collision count > 1
    # missing = collision count == 0
    #

    summary = {
        name: {
            prefix: {
                "wrong": 0,
                "missing": 0,
                "unique": 0,
            }
            for prefix in PREFIXES
        }
        for name, _ in SIGNATURES
    }

    summary["order_pair"] = {
        prefix: {
            "wrong": 0,
            "missing": 0,
            "unique": 0,
        }
        for prefix in PREFIXES
    }

    for target_id, (p, q, n) in enumerate(
        targets,
        1
    ):

        section(
            f"TARGET {target_id:>2} "
            f"({p}, {q})"
        )

        for prefix in PREFIXES:

            # IMPORTANT:
            #
            # Prefix N is evaluated at r = Nth R value.
            #
            # Therefore:
            # prefix 3 -> r=5 -> m=16
            # prefix 4 -> r=7 -> m=22
            # prefix 5 -> r=11 -> m=34
            # prefix 6 -> r=13 -> m=40
            # prefix 7 -> r=17 -> m=52
            #
            r = R_VALUES[prefix - 1]

            m = moduli[r]

            cache = signature_cache[
                (prefix, r)
            ]

            order_table = order_tables[m]

            print(
                f"\nprefix {prefix} "
                f"(r={r}, m={m})"
            )

            for name, _ in SIGNATURES:

                sig = target_sig(
                    p,
                    q,
                    m,
                    name,
                    order_table
                )

                count = cache[name].get(
                    sig,
                    0
                )

                if count == 0:
                    status = "MISSING"
                    summary[name][prefix][
                        "missing"
                    ] += 1

                elif count == 1:
                    status = "UNIQUE"
                    summary[name][prefix][
                        "unique"
                    ] += 1

                else:
                    status = "WRONG"
                    summary[name][prefix][
                        "wrong"
                    ] += 1

                print(
                    f"  {name:<12}: "
                    f"{count:>10,} "
                    f"{status}"
                )

            # Order pair
            name = "order_pair"

            sig = target_sig(
                p,
                q,
                m,
                name,
                order_table
            )

            if sig is None:
                count = 0
            else:
                count = cache[name].get(
                    sig,
                    0
                )

            if count == 0:
                status = "MISSING"
                summary[name][prefix][
                    "missing"
                ] += 1

            elif count == 1:
                status = "UNIQUE"
                summary[name][prefix][
                    "unique"
                ] += 1

            else:
                status = "WRONG"
                summary[name][prefix][
                    "wrong"
                ] += 1

            print(
                f"  {name:<12}: "
                f"{count:>10,} "
                f"{status}"
            )

    # ==========================================================================
    # SUMMARY
    # ==========================================================================

    header(
        "GLOBAL SUMMARY — WRONG COLLISIONS"
    )

    print(
        f"{'signature':<14}"
        + "".join(
            f"{p:>10}"
            for p in PREFIXES
        )
    )

    print("-" * 64)

    all_names = [
        name
        for name, _ in SIGNATURES
    ] + ["order_pair"]

    for name in all_names:

        print(
            f"{name:<14}"
            + "".join(
                f"{summary[name][p]['wrong']:>10}"
                for p in PREFIXES
            )
        )

    # ==========================================================================
    # MISSING SUMMARY
    # ==========================================================================

    header(
        "GLOBAL SUMMARY — TARGET MISSING"
    )

    print(
        f"{'signature':<14}"
        + "".join(
            f"{p:>10}"
            for p in PREFIXES
        )
    )

    print("-" * 64)

    for name in all_names:

        print(
            f"{name:<14}"
            + "".join(
                f"{summary[name][p]['missing']:>10}"
                for p in PREFIXES
            )
        )

    # ==========================================================================
    # UNIQUE SUMMARY
    # ==========================================================================

    header(
        "GLOBAL SUMMARY — ZERO WRONG COLLISIONS"
    )

    print(
        "Count of targets for which the individual "
        "signature has exactly one pair."
    )

    print()

    print(
        f"{'signature':<14}"
        + "".join(
            f"{p:>10}"
            for p in PREFIXES
        )
    )

    print("-" * 64)

    for name in all_names:

        print(
            f"{name:<14}"
            + "".join(
                f"{summary[name][p]['unique']:>10}"
                for p in PREFIXES
            )
        )

    # ==========================================================================
    # FINISH
    # ==========================================================================

    total_time = (
        time.perf_counter()
        - start_total
    )

    header(
        "EXPERIMENT 34 COMPLETE"
    )

    print(
        "This experiment does not perform "
        "pairwise candidate intersections."
    )

    print(
        "It works directly with residue classes, "
        "so the expensive 5,000 choose 2 enumeration "
        "is avoided."
    )

    print()

    print(
        "WRONG    = target signature shared by "
        "multiple prime pairs"
    )

    print(
        "UNIQUE   = exactly one prime pair in pool"
    )

    print(
        "MISSING  = target signature has zero "
        "pairs in the indexed pool"
    )

    print()

    print(
        f"total runtime = "
        f"{total_time:.2f}s"
    )


if __name__ == "__main__":
    main()

