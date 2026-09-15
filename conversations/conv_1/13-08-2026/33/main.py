import random
import time
from collections import defaultdict

import sympy


# ==============================================================================
# KAPPA EXPERIMENT 33
# ULTRA-FAST SIGNATURE DIAGNOSTIC
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
# PRINT HELPERS
# ==============================================================================

def header(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def subheader(title):
    print()
    print("-" * 78)
    print(title)
    print("-" * 78)


# ==============================================================================
# MODULUS CONSTRUCTION
# ==============================================================================

def build_modulus(r):
    """
    Experiment-32 style modulus construction.

    m = 3r + 1
    """

    return 3 * r + 1


def unit_count(m):
    """
    Number of units modulo m.
    """

    return sum(
        1 for x in range(1, m)
        if sympy.gcd(x, m) == 1
    )


# ==============================================================================
# PRIME POOL
# ==============================================================================

def generate_prime_pool():
    primes = []

    while len(primes) < POOL_SIZE:
        x = random.randint(PRIME_MIN, PRIME_MAX)

        if sympy.isprime(x):
            primes.append(x)

    # Remove accidental duplicates while preserving order.
    primes = list(dict.fromkeys(primes))

    # If duplicates reduced the pool, continue.
    while len(primes) < POOL_SIZE:
        x = random.randint(PRIME_MIN, PRIME_MAX)

        if sympy.isprime(x) and x not in primes:
            primes.append(x)

    return primes


# ==============================================================================
# TARGET GENERATION
# ==============================================================================

def generate_targets(prime_pool):
    """
    Generate 12 distinct unordered prime pairs.
    """

    targets = []

    used = set()

    while len(targets) < TARGET_COUNT:
        p, q = random.sample(prime_pool, 2)

        if p > q:
            p, q = q, p

        pair = (p, q)

        if pair in used:
            continue

        used.add(pair)

        n = p * q

        targets.append((p, q, n))

    return targets


# ==============================================================================
# RESIDUE CACHE
# ==============================================================================

def build_residue_cache(primes, moduli):
    cache = {}

    for m in moduli:
        cache[m] = {
            p: p % m
            for p in primes
        }

    return cache


# ==============================================================================
# SIGNATURE FUNCTIONS
# ==============================================================================
#
# These are intentionally kept simple and deterministic.
#
# IMPORTANT:
# The experiment compares candidate pairs by the same functions used for
# the target, so a signature cannot accidentally "lose" the target unless
# the candidate-generation rule itself excludes it.
# ==============================================================================


def Fpair(p, q, m):
    return tuple(sorted((
        p % m,
        q % m,
    )))


def Fsorted(p, q, m):
    return tuple(sorted((
        p % m,
        q % m,
    )))


def Fsum(p, q, m):
    return (p + q) % m


def Fprod(p, q, m):
    return (p * q) % m


def Fdiff(p, q, m):
    return (p - q) % m


def cube_sum(p, q, m):
    return (pow(p, 3, m) + pow(q, 3, m)) % m


def multiplicative_order(a, m):
    """
    Multiplicative order of a modulo m.

    Returns None when gcd(a,m) != 1.
    """

    a %= m

    if sympy.gcd(a, m) != 1:
        return None

    return sympy.n_order(a, m)


def order_pair(p, q, m):
    op = multiplicative_order(p, m)
    oq = multiplicative_order(q, m)

    if op is None or oq is None:
        return None

    return tuple(sorted((op, oq)))


SIGNATURES = [
    ("Fpair", Fpair),
    ("Fsorted", Fsorted),
    ("Fdiff", Fdiff),
    ("Fsum", Fsum),
    ("Fprod", Fprod),
    ("cube_sum", cube_sum),
    ("order_pair", order_pair),
]


# ==============================================================================
# FAST SIGNATURE INDEX
# ==============================================================================

def build_signature_index(prime_pool, m, signature_func, residue_cache):
    """
    Build signature -> count.

    We only need counts for this experiment.
    Candidate pairs are unordered.
    """

    index = defaultdict(int)

    residues = residue_cache[m]

    for i in range(len(prime_pool)):
        p = prime_pool[i]
        rp = residues[p]

        for j in range(i + 1, len(prime_pool)):
            q = prime_pool[j]
            rq = residues[q]

            if signature_func is Fpair:
                sig = tuple(sorted((rp, rq)))

            elif signature_func is Fsorted:
                sig = tuple(sorted((rp, rq)))

            elif signature_func is Fdiff:
                sig = (rp - rq) % m

            elif signature_func is Fsum:
                sig = (rp + rq) % m

            elif signature_func is Fprod:
                sig = (rp * rq) % m

            elif signature_func is cube_sum:
                sig = (
                    pow(rp, 3, m) +
                    pow(rq, 3, m)
                ) % m

            elif signature_func is order_pair:
                op = multiplicative_order(rp, m)
                oq = multiplicative_order(rq, m)

                if op is None or oq is None:
                    continue

                sig = tuple(sorted((op, oq)))

            else:
                sig = signature_func(p, q, m)

            index[sig] += 1

    return index


# ==============================================================================
# TARGET SIGNATURE
# ==============================================================================

def target_signature(p, q, m, signature_func):
    return signature_func(p, q, m)


# ==============================================================================
# DIAGNOSTIC
# ==============================================================================

def run_prefix_diagnostic(
    target_id,
    p,
    q,
    prime_pool,
    m,
    residue_cache,
):
    print()
    print(
        f"prefix {target_id} / m={m} "
        f"target=({p}, {q})"
    )

    prefix_start = time.perf_counter()

    results = {}

    for name, func in SIGNATURES:

        t0 = time.perf_counter()

        index = build_signature_index(
            prime_pool,
            m,
            func,
            residue_cache,
        )

        sig = target_signature(
            p,
            q,
            m,
            func,
        )

        collisions = index.get(sig, 0)

        present = collisions > 0

        elapsed = time.perf_counter() - t0

        results[name] = {
            "collisions": collisions,
            "present": present,
            "time": elapsed,
        }

        status = "YES" if present else "NO "

        print(
            f"  {name:<12}: "
            f"{collisions:>9,} candidates   "
            f"target={status}   "
            f"[{elapsed:.3f}s]"
        )

    total = time.perf_counter() - prefix_start

    print(f"  total prefix time = {total:.3f}s")

    return results


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    experiment_start = time.perf_counter()

    header(
        "KAPPA EXPERIMENT 33\n"
        "ULTRA-FAST SIGNATURE DIAGNOSTIC\n"
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

    # --------------------------------------------------------------------------
    # MODULI
    # --------------------------------------------------------------------------

    subheader("MODULUS INVENTORY")

    moduli = {}

    for r in R_VALUES:
        m = build_modulus(r)
        units = unit_count(m)

        moduli[r] = m

        print(
            f"r={r:>3} "
            f"m={m:>7} "
            f"units={units:>7}"
        )

    # --------------------------------------------------------------------------
    # PRIME POOL
    # --------------------------------------------------------------------------

    subheader("PRIME POOL")

    t0 = time.perf_counter()

    print("Generating prime pool...")

    prime_pool = generate_prime_pool()

    pool_time = time.perf_counter() - t0

    print(f"generated primes = {len(prime_pool):,}")
    print(f"generation time = {pool_time:.3f}s")

    # --------------------------------------------------------------------------
    # TARGETS
    # --------------------------------------------------------------------------

    subheader("TARGETS")

    targets = generate_targets(prime_pool)

    for i, (p, q, n) in enumerate(targets, 1):
        print(
            f"target {i:>2}: "
            f"p={p} "
            f"q={q} "
            f"n={n}"
        )

    # --------------------------------------------------------------------------
    # RESIDUE CACHE
    # --------------------------------------------------------------------------

    subheader("RESIDUE CACHE")

    t0 = time.perf_counter()

    used_moduli = sorted(set(moduli.values()))

    residue_cache = build_residue_cache(
        prime_pool,
        used_moduli,
    )

    print(
        f"residue cache complete "
        f"[{time.perf_counter() - t0:.3f}s]"
    )

    # --------------------------------------------------------------------------
    # DIAGNOSTIC RUN
    # --------------------------------------------------------------------------

    header("SIGNATURE DIAGNOSTIC")

    # Statistics:
    #
    # wrong collision means:
    #   collision count > 1
    #
    # because one candidate is the target itself.
    #
    # missing means:
    #   collision count == 0
    #
    summary = {
        name: {
            prefix: {
                "present": 0,
                "missing": 0,
                "wrong": 0,
            }
            for prefix in PREFIXES
        }
        for name, _ in SIGNATURES
    }

    total_prefixes = 0

    for target_id, (p, q, n) in enumerate(targets, 1):

        subheader(
            f"TARGET {target_id:>2} "
            f"({p}, {q})"
        )

        for prefix in PREFIXES:

            # Prefix means first N R values.
            selected_r = R_VALUES[:prefix]

            # Diagnostic uses the largest selected modulus.
            #
            # This deliberately makes the prefix behavior explicit.
            m = moduli[selected_r[-1]]

            results = run_prefix_diagnostic(
                target_id,
                p,
                q,
                prime_pool,
                m,
                residue_cache,
            )

            total_prefixes += 1

            for name in results:

                count = results[name]["collisions"]

                if count == 0:
                    summary[name][prefix]["missing"] += 1

                else:
                    summary[name][prefix]["present"] += 1

                    if count > 1:
                        summary[name][prefix]["wrong"] += 1

    # --------------------------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------------------------

    header("GLOBAL DIAGNOSTIC SUMMARY")

    print(
        "Each cell = number of targets where the signature had "
        "at least one WRONG collision."
    )

    print()

    print(
        f"{'signature':<14}"
        + "".join(f"{p:>10}" for p in PREFIXES)
    )

    print("-" * 64)

    for name, _ in SIGNATURES:

        row = f"{name:<14}"

        for prefix in PREFIXES:
            row += (
                f"{summary[name][prefix]['wrong']:>10}"
            )

        print(row)

    print()

    print(
        "Each cell = number of targets where the target was "
        "MISSING entirely."
    )

    print()

    print(
        f"{'signature':<14}"
        + "".join(f"{p:>10}" for p in PREFIXES)
    )

    print("-" * 64)

    for name, _ in SIGNATURES:

        row = f"{name:<14}"

        for prefix in PREFIXES:
            row += (
                f"{summary[name][prefix]['missing']:>10}"
            )

        print(row)

    # --------------------------------------------------------------------------
    # IMPORTANT INTERPRETATION
    # --------------------------------------------------------------------------

    header("INTERPRETATION")

    print(
        "This experiment does NOT claim uniqueness."
    )

    print(
        "It only diagnoses whether each signature preserves the "
        "target and how many finite-pool collisions it has."
    )

    print()

    print(
        "A target marked MISSING means the candidate/signature "
        "construction excluded the target itself."
    )

    print(
        "A target with WRONG > 0 has at least one other pair "
        "sharing that signature."
    )

    print(
        "A target with WRONG = 0 has no observed collision for "
        "that individual signature in this finite pool."
    )

    print()

    total_time = time.perf_counter() - experiment_start

    print("=" * 78)
    print("EXPERIMENT 33 COMPLETE")
    print("=" * 78)
    print(f"total runtime = {total_time:.2f}s")
    print()


if __name__ == "__main__":
    main()

