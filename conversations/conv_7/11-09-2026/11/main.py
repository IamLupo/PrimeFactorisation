import math


EXPERIMENT_NUMBER = 33

# Exhaustive range.
# Increase this later if the runtime is acceptable.
MAX_N = 2_000_000


def sieve_primes(limit):
    """Return all primes <= limit."""
    if limit < 2:
        return []

    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0] = 0
    is_prime[1] = 0

    for p in range(2, math.isqrt(limit) + 1):
        if is_prime[p]:
            start = p * p
            count = ((limit - start) // p) + 1

            for multiple in range(start, limit + 1, p):
                is_prime[multiple] = 0

    return [
        n
        for n in range(2, limit + 1)
        if is_prime[n]
    ]


def second_branch_root(s, q):
    """
    Solve:

        4*x + 3*s = 0 (mod q)

    for the least nonnegative x.
    """
    return (-3 * s * pow(4, q - 2, q)) % q


def branch_multiplier(s, q):
    """
    Compute t from:

        4*x2 + 3*s = t*q
    """
    x2 = second_branch_root(s, q)

    value = 4 * x2 + 3 * s

    if value % q != 0:
        raise ValueError("Invalid branch root.")

    return value // q


def factor_semiprime(n, prime_set):
    """
    Return (p,q) only when n is a product of two distinct primes.

    This is used only for exhaustive validation over a small range.
    """
    s = math.isqrt(n)

    for p in prime_set:
        if p > s:
            break

        if n % p != 0:
            continue

        q = n // p

        if q != p and q in prime_set:
            if p < q:
                return p, q
            return q, p

    return None


def construct_semiprimes(max_n, primes):
    """
    Generate every distinct-prime semiprime N <= max_n.

    Returning sorted tuples makes the experiment deterministic.
    """
    prime_set = set(primes)
    cases = []

    for i, p in enumerate(primes):
        if p * p >= max_n:
            break

        for q in primes[i + 1:]:
            n = p * q

            if n > max_n:
                break

            cases.append((n, p, q))

    cases.sort()

    return cases


def fingerprint(n, s, k):
    """
    F_k(N) = (N mod 2^k, s mod 2^k)
    """
    mask = (1 << k) - 1

    return (
        n & mask,
        s & mask,
    )


def analyze_collisions(cases, max_k):
    """
    For each k, find the first fingerprint collision
    with different t.
    """
    results = {}

    for k in range(2, max_k + 1):
        groups = {}
        conflict = None

        for n, p, q in cases:
            s = math.isqrt(n)
            t = branch_multiplier(s, q)

            key = fingerprint(n, s, k)

            previous = groups.get(key)

            if previous is None:
                groups[key] = (n, p, q, t)

            else:
                old_n, old_p, old_q, old_t = previous

                if old_t != t:
                    conflict = {
                        "k": k,
                        "fingerprint": key,
                        "first": (
                            old_n,
                            old_p,
                            old_q,
                            old_t,
                            math.isqrt(old_n),
                        ),
                        "second": (
                            n,
                            p,
                            q,
                            t,
                            s,
                        ),
                    }
                    break

        results[k] = conflict

    return results


def summarize_t_collisions(cases, k):
    """
    Enumerate all conflicting t-pairs for one k.
    """
    groups = {}

    for n, p, q in cases:
        s = math.isqrt(n)
        t = branch_multiplier(s, q)

        key = fingerprint(n, s, k)

        groups.setdefault(key, []).append(
            (n, p, q, t, s)
        )

    pair_counts = {}

    for entries in groups.values():
        t_values = sorted(
            {entry[3] for entry in entries}
        )

        if len(t_values) <= 1:
            continue

        for i in range(len(t_values)):
            for j in range(i + 1, len(t_values)):
                pair = (
                    t_values[i],
                    t_values[j],
                )

                pair_counts[pair] = (
                    pair_counts.get(pair, 0) + 1
                )

    return pair_counts


def run_experiment():
    print(f"START EXPERIMENT {EXPERIMENT_NUMBER}")
    print()

    print("EXHAUSTIVE 2-ADIC COLLISION EXPERIMENT")
    print("--------------------------------------")
    print()
    print("Search every distinct-prime semiprime N <= MAX_N.")
    print()
    print("Fingerprint:")
    print("  F_k(N) = (N mod 2^k, floor(sqrt(N)) mod 2^k)")
    print()
    print("A failure for k means:")
    print("  same fingerprint")
    print("  but different t")
    print()
    print(f"MAX_N = {MAX_N}")
    print()

    primes = sieve_primes(MAX_N)

    print(
        f"Primes generated: {len(primes)}"
    )

    cases = construct_semiprimes(
        MAX_N,
        primes,
    )

    print(
        f"Distinct-prime semiprimes tested: {len(cases)}"
    )

    if not cases:
        print("No cases generated.")
        print()
        print(f"FINISHED EXPERIMENT {EXPERIMENT_NUMBER}")
        return

    print()

    # ---------------------------------------------------------
    # Basic t distribution.
    # ---------------------------------------------------------

    print("ACTUAL t DISTRIBUTION")
    print("---------------------")

    histogram = {}

    for n, p, q in cases:
        s = math.isqrt(n)
        t = branch_multiplier(s, q)

        histogram[t] = histogram.get(t, 0) + 1

    for t in sorted(histogram):
        print(
            f"t={t}: {histogram[t]}"
        )

    print()

    # ---------------------------------------------------------
    # Exhaustive first-conflict search.
    # ---------------------------------------------------------

    max_k = 20

    print("FIRST COUNTEREXAMPLE BY k")
    print("-------------------------")

    collision_results = analyze_collisions(
        cases,
        max_k,
    )

    for k in range(2, max_k + 1):
        conflict = collision_results[k]

        if conflict is None:
            print(
                f"k={k:2d}: "
                f"NO CONFLICT"
            )
            continue

        fp = conflict["fingerprint"]

        n1, p1, q1, t1, s1 = conflict["first"]
        n2, p2, q2, t2, s2 = conflict["second"]

        print(
            f"k={k:2d}: "
            f"CONFLICT "
            f"fingerprint={fp} "
            f"t={t1} vs {t2}"
        )

        print(
            f"       "
            f"N1={n1} "
            f"p1={p1} "
            f"q1={q1} "
            f"s1={s1} "
            f"t1={t1}"
        )

        print(
            f"       "
            f"N2={n2} "
            f"p2={p2} "
            f"q2={q2} "
            f"s2={s2} "
            f"t2={t2}"
        )

    print()

    # ---------------------------------------------------------
    # Stronger summary.
    # ---------------------------------------------------------

    print("NO-CONFLICT SUMMARY")
    print("-------------------")

    for k in range(2, max_k + 1):
        if collision_results[k] is None:
            print(
                f"k={k:2d}: "
                f"no conflicting fingerprint found"
            )

    print()

    # ---------------------------------------------------------
    # Collision-pair structure.
    # ---------------------------------------------------------

    for k in [4, 6, 8, 10, 12]:
        if k > max_k:
            continue

        pair_counts = summarize_t_collisions(
            cases,
            k,
        )

        print(
            f"t COLLISION PAIRS AT k={k}"
        )
        print(
            "-----------------------------"
        )

        if not pair_counts:
            print("No conflicting t-pairs.")
        else:
            for pair in sorted(pair_counts):
                print(
                    f"{pair[0]} vs {pair[1]}: "
                    f"{pair_counts[pair]}"
                )

        print()

    # ---------------------------------------------------------
    # Determine smallest globally collision-free k.
    # ---------------------------------------------------------

    first_clean = None

    for k in range(2, max_k + 1):
        if collision_results[k] is None:
            first_clean = k
            break

    print("FIRST EXHAUSTIVE COLLISION-FREE k")
    print("---------------------------------")

    if first_clean is None:
        print(
            f"No collision-free k <= {max_k}."
        )
    else:
        print(
            f"k={first_clean} is the first "
            f"collision-free value over N <= {MAX_N}."
        )

    print()

    # ---------------------------------------------------------
    # Test the particularly interesting k values.
    # ---------------------------------------------------------

    print("SELECTED k RESULTS")
    print("------------------")

    for k in [4, 6, 8, 10, 12, 14, 16]:
        if k > max_k:
            continue

        conflict = collision_results[k]

        if conflict is None:
            print(
                f"k={k:2d}: "
                f"collision-free"
            )
        else:
            n1, p1, q1, t1, s1 = conflict["first"]
            n2, p2, q2, t2, s2 = conflict["second"]

            print(
                f"k={k:2d}: "
                f"first conflict "
                f"N={n1} (t={t1}) "
                f"vs N={n2} (t={t2})"
            )

    print()

    # ---------------------------------------------------------
    # Look at maximum N represented by each collision-free k.
    # ---------------------------------------------------------

    print("EXHAUSTIVE RANGE")
    print("----------------")

    print(
        f"All semiprimes up to {MAX_N:,} were tested."
    )

    print()

    print("INTERPRETATION")
    print("--------------")
    print("1. This is an exhaustive test, not a random sample.")
    print("2. A collision means identical low bits of N and floor(sqrt(N))")
    print("   but different second-branch multipliers t.")
    print("3. If k=8,10,12 remain collision-free, that is much stronger")
    print("   evidence than the previous random experiments.")
    print("4. Pay particular attention to the first conflicting N values.")
    print("5. The t-pair distribution may reveal which branches are still")
    print("   indistinguishable at low 2-adic precision.")

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT_NUMBER}")


if __name__ == "__main__":
    run_experiment()
