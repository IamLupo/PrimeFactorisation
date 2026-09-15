#!/usr/bin/env python3

import math
import random
from collections import defaultdict

# ============================================================
# KAPPA EXPERIMENT 19
# JOINT SIGNATURE RANKING
#
# No CSV.
# No orbit dumps.
# Only aggregate output.
# ============================================================

SEED = 20260814
random.seed(SEED)

R_VALUES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

TARGETS = 12

# Keep this modest at first.
# The expensive part is the full unit orbit.
MAX_R = 47


# ============================================================
# BASIC NUMBER THEORY
# ============================================================

def gcd(a, b):
    return math.gcd(a, b)


def egcd(a, b):
    if b == 0:
        return (a, 1, 0)

    g, x, y = egcd(b, a % b)
    return (g, y, x - (a // b) * y)


def invmod(a, m):
    g, x, _ = egcd(a, m)
    if g != 1:
        return None
    return x % m


def euler_phi(n):
    result = n
    x = n
    p = 2

    while p * p <= x:
        if x % p == 0:
            while x % p == 0:
                x //= p
            result -= result // p

        p += 1 if p == 2 else 2

    if x > 1:
        result -= result // x

    return result


def multiplicative_order(x, m):
    if gcd(x, m) != 1:
        return 0

    phi = euler_phi(m)
    d = phi

    # factor phi
    factors = []
    z = phi
    p = 2

    while p * p <= z:
        if z % p == 0:
            factors.append(p)
            while z % p == 0:
                z //= p
        p += 1

    if z > 1:
        factors.append(z)

    for p in factors:
        while d % p == 0 and pow(x, d // p, m) == 1:
            d //= p

    return d


# ============================================================
# TARGET GENERATION
# ============================================================

def random_prime(bits):
    while True:
        x = random.getrandbits(bits)
        x |= (1 << (bits - 1))
        x |= 1

        if is_prime(x):
            return x


def is_prime(n):
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

    # Deterministic enough for our ~22-bit primes.
    bases = [2, 3, 5, 7, 11]

    for a in bases:
        if a >= n:
            continue

        x = pow(a, d, n)

        if x in (1, n - 1):
            continue

        for _ in range(s - 1):
            x = pow(x, 2, n)

            if x == n - 1:
                break
        else:
            return False

    return True


def make_target():
    while True:
        p = random_prime(22)
        q = random_prime(22)

        if p == q:
            continue

        n = p * q
        s = p + q
        delta = (q - p) ** 2

        return p, q, n, s, delta


# ============================================================
# MODULUS
# ============================================================

def modulus_for_r(r):
    # This reproduces the second experiment's family:
    #
    # m = r^2 + 2r + 1 = (r+1)^2
    #
    # except r=3 -> 17 etc.?
    #
    # The supplied output corresponds to:
    # m = r^2 + r + 1
    #
    # 2 -> 7
    # 3 -> 13 would be first experiment,
    #
    # so experiment 18 used a different construction.
    #
    # To avoid silently changing your construction, put the
    # exact modulus map below.

    table = {
        2: 7,
        3: 17,
        5: 49,
        7: 97,
        11: 241,
        13: 337,
        17: 577,
        19: 721,
        23: 1057,
        29: 1681,
        31: 1921,
        37: 2737,
        41: 3361,
        43: 3697,
        47: 4417,
    }

    return table[r]


# ============================================================
# ACTUAL FACTOR RESIDUES
# ============================================================

def actual_residues(p, q, m):
    return p % m, q % m


# ============================================================
# UNIT ORBIT
# ============================================================

def units_mod(m):
    return [x for x in range(1, m) if gcd(x, m) == 1]


def orbit_pairs(m, x, y):
    """
    Complete diagonal unit action:

        (x,y) -> (u*x, u*y) mod m

    Returns the complete orbit.
    """

    out = []

    for u in units_mod(m):
        out.append(((u * x) % m, (u * y) % m))

    return out


# ============================================================
# SIGNATURES
# ============================================================

def F(x, m):
    # Cyclotomic-style quadratic signature.
    return (x * x + x + 1) % m


def signature_values(x, y, m):
    fx = F(x, m)
    fy = F(y, m)

    ox = multiplicative_order(x, m)
    oy = multiplicative_order(y, m)

    return {
        "Fpair": (fx, fy),

        "Fsorted": tuple(sorted((fx, fy))),

        "Fsum": (fx + fy) % m,

        "Fprod": (fx * fy) % m,

        "Fdiff": (fx - fy) % m,

        "order_pair": (ox, oy),

        "order_sorted": tuple(sorted((ox, oy))),

        "cube_sum": (pow(x, 3, m) + pow(y, 3, m)) % m,
    }


# ============================================================
# ORBIT FREQUENCY TABLE
# ============================================================

def orbit_frequency_table(m, x, y):
    orbit = orbit_pairs(m, x, y)
    size = len(orbit)

    tables = {
        "Fpair": defaultdict(int),
        "Fsorted": defaultdict(int),
        "Fsum": defaultdict(int),
        "Fprod": defaultdict(int),
        "Fdiff": defaultdict(int),
        "order_pair": defaultdict(int),
        "order_sorted": defaultdict(int),
        "cube_sum": defaultdict(int),
    }

    for a, b in orbit:
        sig = signature_values(a, b, m)

        for name in tables:
            tables[name][sig[name]] += 1

    return tables, size


# ============================================================
# ACTUAL RARITIES
# ============================================================

def actual_rarities(p, q, m):
    x, y = actual_residues(p, q, m)

    tables, orbit_size = orbit_frequency_table(m, x, y)

    actual = signature_values(x, y, m)

    result = {}

    for name, table in tables.items():
        hits = table[actual[name]]
        freq = hits / orbit_size

        result[name] = {
            "hits": hits,
            "orbit": orbit_size,
            "freq": freq,
            "value": actual[name],
        }

    return result


# ============================================================
# JOINT SCORE
# ============================================================

SIGNATURES = [
    "Fpair",
    "Fsorted",
    "Fsum",
    "Fprod",
    "Fdiff",
    "order_pair",
    "order_sorted",
    "cube_sum",
]


def joint_score(all_rarities, signature):
    """
    Product of independent-looking per-modulus rarity.

    Work in log space to avoid underflow.
    """

    log_score = 0.0

    for m in all_rarities:
        f = all_rarities[m][signature]["freq"]

        if f <= 0:
            return float("-inf")

        log_score += math.log10(f)

    return log_score


# ============================================================
# PRINT HELPERS
# ============================================================

def fmt_log(x):
    if x == float("-inf"):
        return "-inf"
    return f"{x:.3f}"


def print_target_header(i, p, q, n, s, delta):
    print()
    print("=" * 78)
    print(f"TARGET {i}")
    print("=" * 78)
    print(f"p       = {p}")
    print(f"q       = {q}")
    print(f"n       = {n}")
    print(f"s       = {s}")
    print(f"Delta   = {delta}")
    print(f"n bits  = {n.bit_length()}")


# ============================================================
# MAIN
# ============================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 19")
    print("JOINT SIGNATURE RANKING")
    print("ACTUAL FACTOR PAIR VS COMPLETE UNIT ORBIT")
    print("=" * 78)

    print(f"random seed = {SEED}")
    print(f"R values    = {R_VALUES}")
    print(f"targets     = {TARGETS}")
    print("CSV output  = NONE")
    print()

    # --------------------------------------------------------
    # MODULUS INVENTORY
    # --------------------------------------------------------

    print("=" * 78)
    print("1. MODULUS INVENTORY")
    print("=" * 78)

    moduli = []

    for r in R_VALUES:
        m = modulus_for_r(r)
        u = euler_phi(m)

        moduli.append(m)

        print(
            f"r={r:>3} "
            f"m={m:>6} "
            f"units={u:>6}"
        )

    # --------------------------------------------------------
    # ACCUMULATORS
    # --------------------------------------------------------

    aggregate_logs = {
        sig: []
        for sig in SIGNATURES
    }

    target_results = []

    # --------------------------------------------------------
    # TARGET LOOP
    # --------------------------------------------------------

    for target_index in range(1, TARGETS + 1):

        p, q, n, s, delta = make_target()

        print_target_header(
            target_index,
            p,
            q,
            n,
            s,
            delta
        )

        all_rarities = {}

        # ----------------------------------------------------
        # Per modulus
        # ----------------------------------------------------

        for m in moduli:

            rarities = actual_rarities(p, q, m)

            all_rarities[m] = rarities

        # ----------------------------------------------------
        # Print compact per-modulus summary
        # ----------------------------------------------------

        print()
        print("PER-MODULUS ACTUAL FREQUENCIES")
        print("-" * 78)

        print(
            f"{'m':>6} "
            f"{'orbit':>7} "
            f"{'Fpair':>10} "
            f"{'Fsum':>10} "
            f"{'Fprod':>10} "
            f"{'Fdiff':>10} "
            f"{'order':>10}"
        )

        print("-" * 78)

        for m in moduli:

            a = all_rarities[m]

            print(
                f"{m:>6} "
                f"{a['Fpair']['orbit']:>7} "
                f"{a['Fpair']['freq']:>10.6g} "
                f"{a['Fsum']['freq']:>10.6g} "
                f"{a['Fprod']['freq']:>10.6g} "
                f"{a['Fdiff']['freq']:>10.6g} "
                f"{a['order_pair']['freq']:>10.6g}"
            )

        # ----------------------------------------------------
        # Joint ranking
        # ----------------------------------------------------

        scores = {}

        for sig in SIGNATURES:
            score = joint_score(all_rarities, sig)
            scores[sig] = score
            aggregate_logs[sig].append(score)

        ranked = sorted(
            scores.items(),
            key=lambda kv: kv[1]
        )

        print()
        print("JOINT RARITY RANKING")
        print("-" * 78)
        print(
            f"{'rank':>4} "
            f"{'signature':<16} "
            f"{'log10(joint freq)':>20}"
        )
        print("-" * 78)

        for rank, (sig, score) in enumerate(ranked, 1):
            print(
                f"{rank:>4} "
                f"{sig:<16} "
                f"{fmt_log(score):>20}"
            )

        # ----------------------------------------------------
        # Best signature
        # ----------------------------------------------------

        best_sig, best_score = ranked[0]

        print()
        print("BEST JOINT SIGNATURE")
        print("-" * 78)
        print(f"signature = {best_sig}")
        print(f"log10 joint frequency = {fmt_log(best_score)}")

        target_results.append({
            "best": best_sig,
            "score": best_score,
        })

    # ========================================================
    # FINAL SUMMARY
    # ========================================================

    print()
    print("=" * 78)
    print("FINAL SUMMARY")
    print("=" * 78)

    print()
    print("Median log10(joint frequency)")
    print("-" * 78)

    medians = {}

    for sig in SIGNATURES:
        values = sorted(aggregate_logs[sig])

        k = len(values)

        if k % 2:
            med = values[k // 2]
        else:
            med = (
                values[k // 2 - 1]
                + values[k // 2]
            ) / 2

        medians[sig] = med

    for sig, value in sorted(
        medians.items(),
        key=lambda kv: kv[1]
    ):
        print(
            f"{sig:<16} "
            f"{fmt_log(value):>12}"
        )

    # --------------------------------------------------------
    # Win counts
    # --------------------------------------------------------

    print()
    print("JOINT-RANK WIN COUNTS")
    print("-" * 78)

    wins = defaultdict(int)

    for row in target_results:
        wins[row["best"]] += 1

    for sig in SIGNATURES:
        print(
            f"{sig:<16} "
            f"{wins[sig]:>4}/{TARGETS}"
        )

    # --------------------------------------------------------
    # Most important comparison
    # --------------------------------------------------------

    print()
    print("KEY COMPARISON")
    print("-" * 78)

    fpair_med = medians["Fpair"]

    for sig in [
        "Fsorted",
        "Fsum",
        "Fprod",
        "Fdiff",
        "order_pair",
        "order_sorted",
        "cube_sum",
    ]:

        improvement = medians[sig] - fpair_med

        print(
            f"{sig:<16} "
            f"relative log10 score vs Fpair = "
            f"{improvement:+.3f}"
        )

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()

