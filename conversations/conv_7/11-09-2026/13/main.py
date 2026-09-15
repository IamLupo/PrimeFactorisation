import math
import random


EXPERIMENT_NUMBER = 35


def is_probable_prime(n):
    """Deterministic Miller-Rabin for 64-bit integers."""
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

        composite = True

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                composite = False
                break

        if composite:
            return False

    return True


def random_prime(bits, rng):
    """Generate a random prime with exactly 'bits' bits."""
    while True:
        candidate = rng.getrandbits(bits)

        candidate |= 1 << (bits - 1)
        candidate |= 1

        if is_probable_prime(candidate):
            return candidate


def second_branch_root(s, q):
    """Solve 4*x2 + 3*s = 0 (mod q)."""
    return (-3 * s * pow(4, q - 2, q)) % q


def branch_multiplier(s, q, x2):
    """Compute t from 4*x2 + 3*s = t*q."""
    value = 4 * x2 + 3 * s

    if value % q != 0:
        raise ValueError("Invalid second-branch root.")

    return value // q


def construct_case(p, q):
    """Construct one factorization case."""
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
        "x2": x2,
        "t": t,
    }


def generate_cases(bits_p, bits_q, count, rng):
    """Generate distinct semiprimes."""
    cases = []
    seen = set()

    while len(cases) < count:
        p = random_prime(bits_p, rng)
        q = random_prime(bits_q, rng)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n in seen:
            continue

        seen.add(n)
        cases.append(construct_case(p, q))

    return cases


def fingerprint(case, k):
    """Return (N mod 2^k, s mod 2^k)."""
    mask = (1 << k) - 1
    return (
        case["N"] & mask,
        case["s"] & mask,
    )


def v2(n):
    """
    2-adic valuation of nonzero integer n.

    v2(0) returns infinity.
    """
    n = abs(n)

    if n == 0:
        return float("inf")

    return (n & -n).bit_length() - 1


def find_35_collisions(cases, k):
    """
    Find fingerprint collisions where one case has t=3
    and the other has t=5.
    """
    groups = {}

    for case in cases:
        if case["t"] not in (3, 5):
            continue

        key = fingerprint(case, k)
        groups.setdefault(key, []).append(case)

    collisions = []

    for key, entries in groups.items():
        has3 = [c for c in entries if c["t"] == 3]
        has5 = [c for c in entries if c["t"] == 5]

        if not has3 or not has5:
            continue

        for c3 in has3:
            for c5 in has5:
                collisions.append(
                    (key, c3, c5)
                )

    return collisions


def analyze_collision(c3, c5):
    """Calculate exact delta and 2-adic information."""
    delta_s = c5["s"] - c3["s"]
    delta_n = c5["N"] - c3["N"]
    delta_d = c5["D"] - c3["D"]
    delta_q = c5["q"] - c3["q"]

    # The branch coordinates satisfy:
    #
    # t=3: x3 = 3(q-s)/4
    # t=5: x5 = (5q-3s)/4
    #
    # The same q is not assumed between collision cases,
    # so we measure both individually and their differences.

    x3 = c3["x2"]
    x5 = c5["x2"]

    return {
        "delta_s": delta_s,
        "delta_n": delta_n,
        "delta_d": delta_d,
        "delta_q": delta_q,
        "v2_delta_s": v2(delta_s),
        "v2_delta_n": v2(delta_n),
        "v2_delta_d": v2(delta_d),
        "v2_delta_q": v2(delta_q),
        "x3": x3,
        "x5": x5,
        "x5_minus_x3": x5 - x3,
        "q3_half": c3["q"] // 2,
        "q5_half": c5["q"] // 2,
    }


def run_experiment():
    print(f"START EXPERIMENT {EXPERIMENT_NUMBER}")
    print()

    print("t=3 vs t=5 COLLISION EXPERIMENT")
    print("--------------------------------")
    print()
    print("Study collisions where:")
    print("  F_k(N1) = F_k(N2)")
    print("  t(N1)=3")
    print("  t(N2)=5")
    print()
    print("For each collision calculate:")
    print("  delta N")
    print("  delta s")
    print("  delta D")
    print("  delta q")
    print("  v2(delta N)")
    print("  v2(delta s)")
    print("  v2(delta D)")
    print("  v2(delta q)")
    print()
    print("Goal:")
    print("  Find an exact 2-adic structure behind the")
    print("  surviving 3 <-> 5 ambiguity.")
    print()

    rng = random.Random(35)

    # Large enough that s is well above 4096.
    scales = [
        (18, 18, 3000),
        (22, 22, 3000),
        (26, 26, 3000),
        (30, 30, 3000),
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

        print(
            f"bits={bits_p}+{bits_q} "
            f"cases={len(cases)} "
            f"s_min={min(c['s'] for c in cases)} "
            f"s_max={max(c['s'] for c in cases)}"
        )

    print()
    print(f"Total cases: {len(all_cases)}")
    print()

    # ---------------------------------------------------------
    # Find collisions at different k.
    # ---------------------------------------------------------

    print("3 <-> 5 COLLISION COUNTS")
    print("-------------------------")

    collision_sets = {}

    for k in [6, 8, 10, 12, 14, 16, 18]:
        collisions = find_35_collisions(
            all_cases,
            k,
        )

        collision_sets[k] = collisions

        print(
            f"k={k:2d}: "
            f"{len(collisions)} collision pairs"
        )

    print()

    # ---------------------------------------------------------
    # Aggregate v2 statistics.
    # ---------------------------------------------------------

    print("2-ADIC DELTA STATISTICS")
    print("-----------------------")

    for k in [8, 10, 12]:
        collisions = collision_sets[k]

        if not collisions:
            print(f"k={k}: no collisions")
            continue

        v2_s = []
        v2_n = []
        v2_d = []
        v2_q = []

        for _, c3, c5 in collisions:
            a = analyze_collision(c3, c5)

            v2_s.append(a["v2_delta_s"])
            v2_n.append(a["v2_delta_n"])
            v2_d.append(a["v2_delta_d"])
            v2_q.append(a["v2_delta_q"])

        print(
            f"k={k}: "
            f"count={len(collisions)} "
            f"min_v2(dS)={min(v2_s)} "
            f"min_v2(dN)={min(v2_n)} "
            f"min_v2(dD)={min(v2_d)} "
            f"min_v2(dQ)={min(v2_q)}"
        )

    print()

    # ---------------------------------------------------------
    # Check whether observed delta_s is always divisible by 2^k.
    # ---------------------------------------------------------

    print("FINGERPRINT CONSEQUENCES")
    print("------------------------")

    for k in [8, 10, 12]:
        collisions = collision_sets[k]

        if not collisions:
            continue

        all_s_divisible = all(
            a["s"] % (1 << k) == b["s"] % (1 << k)
            for _, a, b in collisions
        )

        all_n_divisible = all(
            a["N"] % (1 << k) == b["N"] % (1 << k)
            for _, a, b in collisions
        )

        print(
            f"k={k}: "
            f"same_s_mod_2^k={all_s_divisible} "
            f"same_N_mod_2^k={all_n_divisible}"
        )

    print()

    # ---------------------------------------------------------
    # Search for exact relationships between deltas.
    # ---------------------------------------------------------

    print("DELTA RELATION TESTS")
    print("--------------------")

    for k in [8, 10, 12]:
        collisions = collision_sets[k]

        if not collisions:
            continue

        tests = {
            "dN-ds": 0,
            "dD-ds": 0,
            "dD-dN": 0,
        }

        for _, c3, c5 in collisions:
            a = analyze_collision(c3, c5)

            if a["delta_n"] != a["delta_s"]:
                tests["dN-ds"] += 1

            if a["delta_d"] != a["delta_s"]:
                tests["dD-ds"] += 1

            if a["delta_d"] != a["delta_n"]:
                tests["dD-dN"] += 1

        print(
            f"k={k}: "
            f"dN!=ds={tests['dN-ds']}/{len(collisions)} "
            f"dD!=ds={tests['dD-ds']}/{len(collisions)} "
            f"dD!=dN={tests['dD-dN']}/{len(collisions)}"
        )

    print()

    # ---------------------------------------------------------
    # Test the branch-coordinate identities.
    # ---------------------------------------------------------

    print("BRANCH COORDINATE IDENTITIES")
    print("----------------------------")

    identity_failures = 0
    identity_total = 0

    for _, c3, c5 in collision_sets.get(12, []):
        identity_total += 1

        # For t=3:
        #   4*x3 = 3(q3-s3)
        lhs3 = 4 * c3["x2"]
        rhs3 = 3 * (
            c3["q"] - c3["s"]
        )

        # For t=5:
        #   4*x5 = 5q5-3s5
        lhs5 = 4 * c5["x2"]
        rhs5 = 5 * c5["q"] - 3 * c5["s"]

        if lhs3 != rhs3 or lhs5 != rhs5:
            identity_failures += 1

    print(
        f"k=12 branch identities: "
        f"{identity_total - identity_failures}/"
        f"{identity_total}"
    )

    print()

    # ---------------------------------------------------------
    # Detailed k=12 collisions.
    # ---------------------------------------------------------

    print("DETAILED k=12 COLLISIONS")
    print("------------------------")

    collisions = collision_sets[12]

    if not collisions:
        print("No k=12 collisions found.")

    for index, (fp, c3, c5) in enumerate(
        collisions[:30],
        start=1,
    ):
        info = analyze_collision(
            c3,
            c5,
        )

        print(
            f"Collision {index}"
        )

        print(
            f"  fingerprint={fp}"
        )

        print(
            f"  t=3: "
            f"N={c3['N']} "
            f"p={c3['p']} "
            f"q={c3['q']} "
            f"s={c3['s']} "
            f"D={c3['D']} "
            f"x2={c3['x2']}"
        )

        print(
            f"  t=5: "
            f"N={c5['N']} "
            f"p={c5['p']} "
            f"q={c5['q']} "
            f"s={c5['s']} "
            f"D={c5['D']} "
            f"x2={c5['x2']}"
        )

        print(
            f"  deltaN={info['delta_n']}"
        )

        print(
            f"  deltas={info['delta_s']}"
        )

        print(
            f"  deltaD={info['delta_d']}"
        )

        print(
            f"  deltaq={info['delta_q']}"
        )

        print(
            f"  v2(deltaN)={info['v2_delta_n']}"
        )

        print(
            f"  v2(deltas)={info['v2_delta_s']}"
        )

        print(
            f"  v2(deltaD)={info['v2_delta_d']}"
        )

        print(
            f"  v2(deltaq)={info['v2_delta_q']}"
        )

        print()

    # ---------------------------------------------------------
    # Compare x2 values.
    # ---------------------------------------------------------

    print("x2 COMPARISON")
    print("-------------")

    for k in [8, 10, 12]:
        collisions = collision_sets[k]

        if not collisions:
            continue

        differences = [
            abs(
                c5["x2"] - c3["x2"]
            )
            for _, c3, c5 in collisions
        ]

        print(
            f"k={k}: "
            f"min|x5-x3|={min(differences)} "
            f"max|x5-x3|={max(differences)}"
        )

    print()

    # ---------------------------------------------------------
    # Find repeated normalized delta patterns.
    # ---------------------------------------------------------

    print("REPEATED v2(delta) PATTERNS")
    print("---------------------------")

    for k in [8, 10, 12]:
        collisions = collision_sets[k]

        patterns = {}

        for _, c3, c5 in collisions:
            a = analyze_collision(
                c3,
                c5,
            )

            pattern = (
                a["v2_delta_n"],
                a["v2_delta_s"],
                a["v2_delta_d"],
                a["v2_delta_q"],
            )

            patterns[pattern] = (
                patterns.get(pattern, 0) + 1
            )

        print(
            f"k={k}: "
            f"distinct patterns={len(patterns)}"
        )

        for pattern, count in sorted(
            patterns.items(),
            key=lambda item: item[1],
            reverse=True,
        )[:10]:
            print(
                f"  {pattern}: {count}"
            )

    print()

    print("INTERPRETATION")
    print("--------------")
    print("1. Focus on the k=12 collisions.")
    print("2. Check whether delta_s, delta_N and delta_D")
    print("   have a common higher 2-adic structure.")
    print("3. Check whether the 3<->5 ambiguity produces")
    print("   a characteristic delta pattern.")
    print("4. Compare x2_3 and x2_5 as well.")
    print("5. The goal is to derive an exact equation for")
    print("   why the same low-bit fingerprint can support")
    print("   both t=3 and t=5.")

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT_NUMBER}")


if __name__ == "__main__":
    run_experiment()
