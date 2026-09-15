import math
import random


EXPERIMENT_NUMBER = 34


def is_probable_prime(n):
    """
    Deterministic Miller-Rabin for 64-bit integers.
    """
    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small_primes:
        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    bases = (
        2,
        325,
        9375,
        28178,
        450775,
        9780504,
        1795265022,
    )

    for a in bases:
        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        witness = True

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                witness = False
                break

        if witness:
            return False

    return True


def random_prime(bits, rng):
    """Generate a random prime with exactly 'bits' bits."""
    if bits < 3:
        raise ValueError("bits must be >= 3")

    while True:
        candidate = rng.getrandbits(bits)

        candidate |= 1 << (bits - 1)
        candidate |= 1

        if is_probable_prime(candidate):
            return candidate


def generate_semiprime(bits_p, bits_q, rng):
    """Generate distinct primes p < q."""
    while True:
        p = random_prime(bits_p, rng)
        q = random_prime(bits_q, rng)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        return p, q


def second_branch_root(s, q):
    """
    Solve:

        4*x2 + 3*s = 0 (mod q)
    """
    return (-3 * s * pow(4, q - 2, q)) % q


def branch_multiplier(s, q, x2):
    """Compute t from 4*x2 + 3*s = t*q."""
    value = 4 * x2 + 3 * s

    if value % q != 0:
        raise ValueError("Invalid second-branch root.")

    return value // q


def construct_case(p, q):
    """Construct all quantities needed by the experiment."""
    n = p * q
    s = math.isqrt(n)
    d = n - s * s

    x2 = second_branch_root(s, q)
    t = branch_multiplier(s, q, x2)

    return {
        "N": n,
        "p": p,
        "q": q,
        "s": s,
        "D": d,
        "t": t,
    }


def generate_cases(bits_p, bits_q, count, rng):
    """Generate distinct semiprimes."""
    cases = []
    seen = set()

    while len(cases) < count:
        p, q = generate_semiprime(
            bits_p,
            bits_q,
            rng,
        )

        n = p * q

        if n in seen:
            continue

        seen.add(n)

        case = construct_case(p, q)

        if case["s"] <= 4096:
            continue

        cases.append(case)

    return cases


def fingerprint(case, k):
    """
    F_k(N) = (N mod 2^k, s mod 2^k)
    """
    mask = (1 << k) - 1

    return (
        case["N"] & mask,
        case["s"] & mask,
    )


def find_first_collision(cases, k):
    """
    Find first collision with different t.
    """
    groups = {}

    for case in cases:
        key = fingerprint(case, k)

        previous = groups.get(key)

        if previous is None:
            groups[key] = case
            continue

        if previous["t"] != case["t"]:
            return previous, case, key

    return None


def find_all_conflicting_pairs(cases, k):
    """
    Return all fingerprint groups containing multiple t values.
    """
    groups = {}

    for case in cases:
        key = fingerprint(case, k)

        groups.setdefault(key, []).append(case)

    conflicts = []

    for key, entries in groups.items():
        t_values = sorted({
            entry["t"]
            for entry in entries
        })

        if len(t_values) <= 1:
            continue

        conflicts.append(
            (
                key,
                entries,
                t_values,
            )
        )

    return conflicts


def print_collision(label, first, second, fingerprint_value, k):
    """Print detailed collision information."""
    n1 = first["N"]
    n2 = second["N"]

    s1 = first["s"]
    s2 = second["s"]

    d1 = first["D"]
    d2 = second["D"]

    print(
        f"{label} k={k}"
    )

    print(
        f"  fingerprint={fingerprint_value}"
    )

    print(
        f"  N1={n1} "
        f"p1={first['p']} "
        f"q1={first['q']} "
        f"s1={s1} "
        f"D1={d1} "
        f"t1={first['t']}"
    )

    print(
        f"  N2={n2} "
        f"p2={second['p']} "
        f"q2={second['q']} "
        f"s2={s2} "
        f"D2={d2} "
        f"t2={second['t']}"
    )

    print(
        f"  deltaN={n2 - n1}"
    )

    print(
        f"  deltas={s2 - s1}"
    )

    print(
        f"  deltaD={d2 - d1}"
    )

    print(
        f"  same_s={s1 == s2}"
    )


def run_experiment():
    print(f"START EXPERIMENT {EXPERIMENT_NUMBER}")
    print()

    print("LARGE-s 2-ADIC COLLISION EXPERIMENT")
    print("-----------------------------------")
    print()
    print("Previous exhaustive test had s < 4096.")
    print()
    print("This experiment deliberately enforces:")
    print("  s = floor(sqrt(N)) > 4096")
    print()
    print("Fingerprint:")
    print("  F_k(N) = (N mod 2^k, s mod 2^k)")
    print()
    print("Goal:")
    print("  Search for the first k=12 collision after s")
    print("  itself is no longer represented exactly.")
    print()

    rng = random.Random(34)

    # Several scales. All produce s >> 4096.
    scales = [
        (18, 18, 2500),
        (22, 22, 2500),
        (26, 26, 2500),
        (30, 30, 2500),
    ]

    all_cases = []

    print("DATASET SCALES")
    print("--------------")

    for bits_p, bits_q, count in scales:
        cases = generate_cases(
            bits_p,
            bits_q,
            count,
            rng,
        )

        all_cases.extend(cases)

        s_values = [
            case["s"]
            for case in cases
        ]

        print(
            f"p,q bits={bits_p},{bits_q} "
            f"cases={len(cases)} "
            f"s_min={min(s_values)} "
            f"s_max={max(s_values)}"
        )

    print()
    print(
        f"Total cases: {len(all_cases)}"
    )

    print()

    # ---------------------------------------------------------
    # Basic t distribution.
    # ---------------------------------------------------------

    print("ACTUAL t DISTRIBUTION")
    print("---------------------")

    histogram = {}

    for case in all_cases:
        t = case["t"]
        histogram[t] = histogram.get(t, 0) + 1

    for t in sorted(histogram):
        print(
            f"t={t}: {histogram[t]}"
        )

    print()

    # ---------------------------------------------------------
    # Collision test.
    # ---------------------------------------------------------

    print("FIRST COLLISION BY k")
    print("--------------------")

    first_collisions = {}

    for k in [8, 10, 12, 14, 16]:
        collision = find_first_collision(
            all_cases,
            k,
        )

        first_collisions[k] = collision

        if collision is None:
            print(
                f"k={k:2d}: NO CONFLICT"
            )
        else:
            first, second, fp = collision

            print(
                f"k={k:2d}: "
                f"CONFLICT "
                f"fingerprint={fp} "
                f"t={first['t']} vs {second['t']}"
            )

    print()

    # ---------------------------------------------------------
    # Detailed collisions.
    # ---------------------------------------------------------

    for k in [8, 10, 12, 14, 16]:
        collision = first_collisions[k]

        if collision is None:
            continue

        first, second, fp = collision

        print_collision(
            "FIRST CONFLICT",
            first,
            second,
            fp,
            k,
        )

        print()

    # ---------------------------------------------------------
    # Collision counts.
    # ---------------------------------------------------------

    print("CONFLICT COUNTS")
    print("---------------")

    for k in [8, 10, 12, 14, 16]:
        conflicts = find_all_conflicting_pairs(
            all_cases,
            k,
        )

        print(
            f"k={k:2d}: "
            f"conflicting fingerprint groups="
            f"{len(conflicts)}"
        )

    print()

    # ---------------------------------------------------------
    # t-pair structure.
    # ---------------------------------------------------------

    print("t COLLISION PAIRS")
    print("-----------------")

    for k in [8, 10, 12, 14, 16]:
        conflicts = find_all_conflicting_pairs(
            all_cases,
            k,
        )

        pair_counts = {}

        for _, entries, t_values in conflicts:
            for i in range(len(t_values)):
                for j in range(i + 1, len(t_values)):
                    pair = (
                        t_values[i],
                        t_values[j],
                    )

                    pair_counts[pair] = (
                        pair_counts.get(pair, 0) + 1
                    )

        print(
            f"k={k}:"
        )

        if not pair_counts:
            print(
                "  no conflicting t-pairs"
            )
        else:
            for pair in sorted(pair_counts):
                print(
                    f"  {pair[0]} vs {pair[1]}: "
                    f"{pair_counts[pair]}"
                )

    print()

    # ---------------------------------------------------------
    # Same-s collisions.
    #
    # If s1=s2, then identical N mod 2^k means
    # identical D mod 2^k because:
    #
    # N=s^2+D.
    # ---------------------------------------------------------

    print("SAME-s COLLISION ANALYSIS")
    print("-------------------------")

    for k in [8, 10, 12, 14, 16]:
        conflicts = find_all_conflicting_pairs(
            all_cases,
            k,
        )

        same_s = 0
        different_s = 0

        for _, entries, _ in conflicts:
            for i in range(len(entries)):
                for j in range(i + 1, len(entries)):
                    if entries[i]["t"] == entries[j]["t"]:
                        continue

                    if entries[i]["s"] == entries[j]["s"]:
                        same_s += 1
                    else:
                        different_s += 1

        print(
            f"k={k:2d}: "
            f"same-s={same_s} "
            f"different-s={different_s}"
        )

    print()

    # ---------------------------------------------------------
    # Analyze delta_s.
    # ---------------------------------------------------------

    print("DELTA-s STRUCTURE OF CONFLICTS")
    print("------------------------------")

    for k in [8, 10, 12]:
        conflicts = find_all_conflicting_pairs(
            all_cases,
            k,
        )

        deltas = []

        for _, entries, _ in conflicts:
            if len(entries) < 2:
                continue

            # Only compare different-t entries.
            for i in range(len(entries)):
                for j in range(i + 1, len(entries)):
                    if entries[i]["t"] == entries[j]["t"]:
                        continue

                    delta_s = abs(
                        entries[i]["s"]
                        - entries[j]["s"]
                    )

                    deltas.append(delta_s)

        if not deltas:
            print(
                f"k={k}: no conflicts"
            )
            continue

        powers = [
            value
            for value in deltas
            if value > 0
            and (
                value
                & (value - 1)
            ) == 0
        ]

        print(
            f"k={k}: "
            f"conflicts={len(deltas)} "
            f"min|delta_s|={min(deltas)} "
            f"max|delta_s|={max(deltas)} "
            f"power-of-two delta_s={len(powers)}"
        )

    print()

    # ---------------------------------------------------------
    # Test larger k on the entire dataset.
    # ---------------------------------------------------------

    print("LARGER-k STRESS TEST")
    print("--------------------")

    for k in range(12, 21):
        collision = find_first_collision(
            all_cases,
            k,
        )

        if collision is None:
            print(
                f"k={k:2d}: "
                f"collision-free"
            )
        else:
            first, second, fp = collision

            print(
                f"k={k:2d}: "
                f"conflict "
                f"N1={first['N']} "
                f"N2={second['N']} "
                f"t={first['t']} vs {second['t']}"
            )

    print()

    # ---------------------------------------------------------
    # Explicitly locate the first k=12 conflict if one exists.
    # ---------------------------------------------------------

    print("k=12 VERDICT")
    print("-------------")

    collision = first_collisions[12]

    if collision is None:
        print(
            "No k=12 collision found in the "
            "large-s dataset."
        )
    else:
        first, second, fp = collision

        print(
            "k=12 DOES NOT remain collision-free."
        )

        print_collision(
            "COUNTEREXAMPLE",
            first,
            second,
            fp,
            12,
        )

    print()

    print("INTERPRETATION")
    print("--------------")
    print("1. This dataset deliberately has s > 4096.")
    print("2. Therefore s mod 4096 is no longer the exact s.")
    print("3. A k=12 collision here is much stronger evidence")
    print("   against a fixed 12-bit fingerprint.")
    print("4. If k=12 remains collision-free, that is substantially")
    print("   more interesting than the previous exhaustive result.")
    print("5. Examine same-s versus different-s collisions and")
    print("   the delta_s/delta_D structure.")

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT_NUMBER}")


if __name__ == "__main__":
    run_experiment()
