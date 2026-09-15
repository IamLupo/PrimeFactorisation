import math
import random


EXPERIMENT_NUMBER = 30


def sieve_primes(limit):
    """Return all primes <= limit."""
    if limit < 2:
        return []

    is_prime = [True] * (limit + 1)
    is_prime[0] = False
    is_prime[1] = False

    for p in range(2, math.isqrt(limit) + 1):
        if not is_prime[p]:
            continue

        for multiple in range(p * p, limit + 1, p):
            is_prime[multiple] = False

    return [
        n
        for n in range(2, limit + 1)
        if is_prime[n]
    ]


def second_branch_root(s, q):
    """
    Solve:

        4*x + 3*s = 0 (mod q)

    using the least nonnegative residue.
    """
    return (-3 * s * pow(4, q - 2, q)) % q


def actual_t(s, q, x):
    """Compute t from 4*x + 3*s = t*q."""
    value = 4 * x + 3 * s

    if value % q != 0:
        raise ValueError(
            "Invalid second-branch root: "
            f"4*{x} + 3*{s} is not divisible by {q}"
        )

    return value // q


def full_t_candidates(s, q):
    """
    Return t values satisfying the complete conditions:

        0 <= x2 < q

    with

        x2 = (t*q - 3*s)/4.

    Equivalently:

        3*s <= t*q < 4*q + 3*s
        t*q == 3*s (mod 4)
    """
    candidates = []

    for t in range(1, 7):
        value = t * q

        if value < 3 * s:
            continue

        if value >= 4 * q + 3 * s:
            continue

        if (value - 3 * s) % 4 != 0:
            continue

        candidates.append(t)

    return candidates


def weak_t_candidates(s, q):
    """
    The weaker test used earlier:

        t*q >= 3*s
        t*q == 3*s (mod 4)

    without the x2 < q upper bound.
    """
    candidates = []

    for t in range(1, 7):
        value = t * q

        if value < 3 * s:
            continue

        if (value - 3 * s) % 4 != 0:
            continue

        candidates.append(t)

    return candidates


def generate_semiprimes(primes, count, seed):
    """Generate distinct semiprime factor pairs."""
    random.seed(seed)

    pairs = set()

    while len(pairs) < count:
        p = random.choice(primes[10:700])
        q = random.choice(primes[10:700])

        if p == q:
            continue

        if p > q:
            p, q = q, p

        pairs.add((p, q))

    return sorted(pairs)


def analyze_case(p, q):
    """Analyze one semiprime."""
    n = p * q
    s = math.isqrt(n)
    d = n - s * s

    x2 = second_branch_root(s, q)
    t = actual_t(s, q, x2)

    weak_candidates = weak_t_candidates(s, q)
    full_candidates = full_t_candidates(s, q)

    xp = s + 1 - p
    y = p + q - 2 * s - 1

    return {
        "N": n,
        "p": p,
        "q": q,
        "s": s,
        "D": d,
        "x2": x2,
        "t": t,
        "weak_candidates": weak_candidates,
        "full_candidates": full_candidates,
        "xp": xp,
        "y": y,
        "N_mod4": n % 4,
        "s_mod4": s % 4,
        "D_mod4": d % 4,
        "q_mod4": q % 4,
    }


def run_experiment():
    print(f"START EXPERIMENT {EXPERIMENT_NUMBER}")
    print()

    print("REVERSE t-BRANCH EXPERIMENT")
    print("---------------------------")
    print()
    print("Second branch:")
    print("  4*x2 + 3*s = t*q")
    print()
    print("Full conditions:")
    print("  3*s <= t*q")
    print("  t*q < 4*q + 3*s")
    print("  t*q == 3*s (mod 4)")
    print()
    print("Goal:")
    print("  Determine how much of the correct t-branch")
    print("  can be identified from arithmetic information.")
    print()

    primes = sieve_primes(5000)

    cases = generate_semiprimes(
        primes,
        1000,
        seed=30
    )

    analyses = [
        analyze_case(p, q)
        for p, q in cases
    ]

    total = len(analyses)

    print("AGGREGATE RESULTS")
    print("-----------------")

    actual_in_full = sum(
        a["t"] in a["full_candidates"]
        for a in analyses
    )

    unique_full = sum(
        len(a["full_candidates"]) == 1
        for a in analyses
    )

    print(
        f"Actual t in full candidates      : "
        f"{actual_in_full}/{total}"
    )

    print(
        f"Exactly one full candidate        : "
        f"{unique_full}/{total}"
    )

    print()

    print("WEAK VS FULL CANDIDATE SETS")
    print("---------------------------")

    weak_two = 0
    weak_one = 0
    weak_more = 0

    for a in analyses:
        count = len(a["weak_candidates"])

        if count == 1:
            weak_one += 1
        elif count == 2:
            weak_two += 1
        else:
            weak_more += 1

    print(f"Weak: 1 candidate                 : {weak_one}/{total}")
    print(f"Weak: 2 candidates                : {weak_two}/{total}")
    print(f"Weak: >2 candidates               : {weak_more}/{total}")

    print()

    # ---------------------------------------------------------
    # t distribution
    # ---------------------------------------------------------

    print("ACTUAL t DISTRIBUTION")
    print("---------------------")

    histogram = {}

    for a in analyses:
        t = a["t"]
        histogram[t] = histogram.get(t, 0) + 1

    for t in sorted(histogram):
        print(
            f"t={t}: "
            f"{histogram[t]:4d}/{total}"
        )

    print()

    # ---------------------------------------------------------
    # Compare weak and full candidates
    # ---------------------------------------------------------

    print("FULL-BOUND DISAMBIGUATION")
    print("-------------------------")

    disambiguated = 0

    for a in analyses:
        if (
            len(a["weak_candidates"]) > 1
            and len(a["full_candidates"]) == 1
        ):
            disambiguated += 1

    print(
        f"Two weak candidates -> one full candidate: "
        f"{disambiguated}/{total}"
    )

    print()

    # ---------------------------------------------------------
    # Candidate counts
    # ---------------------------------------------------------

    print("FULL CANDIDATE COUNTS")
    print("---------------------")

    full_histogram = {}

    for a in analyses:
        count = len(a["full_candidates"])
        full_histogram[count] = (
            full_histogram.get(count, 0) + 1
        )

    for count in sorted(full_histogram):
        print(
            f"{count} candidate(s): "
            f"{full_histogram[count]:4d}/{total}"
        )

    print()

    # ---------------------------------------------------------
    # Residue groups
    # ---------------------------------------------------------

    print("N MOD 4 / s MOD 4 / D MOD 4 GROUPS")
    print("-----------------------------------")

    groups = {}

    for a in analyses:
        key = (
            a["N_mod4"],
            a["s_mod4"],
            a["D_mod4"],
        )

        groups.setdefault(key, []).append(a)

    for key in sorted(groups):
        subset = groups[key]

        t_values = {}
        q_values = {}

        for a in subset:
            t_values[a["t"]] = (
                t_values.get(a["t"], 0) + 1
            )

            q_values[a["q_mod4"]] = (
                q_values.get(a["q_mod4"], 0) + 1
            )

        t_distribution = " ".join(
            f"t{t}={t_values[t]}"
            for t in sorted(t_values)
        )

        q_distribution = " ".join(
            f"q4={r}:{q_values[r]}"
            for r in sorted(q_values)
        )

        print(
            f"Nmod4={key[0]} "
            f"smod4={key[1]} "
            f"Dmod4={key[2]} "
            f"cases={len(subset)} "
            f"{t_distribution} "
            f"[{q_distribution}]"
        )

    print()

    # ---------------------------------------------------------
    # Examples where weak test gives two candidates
    # ---------------------------------------------------------

    print("WEAKLY AMBIGUOUS CASES")
    print("-----------------------")

    shown = 0

    for a in analyses:
        if len(a["weak_candidates"]) != 2:
            continue

        print(
            f"N={a['N']} "
            f"p={a['p']} "
            f"q={a['q']} "
            f"s={a['s']} "
            f"x2={a['x2']} "
            f"actual_t={a['t']} "
            f"weak={a['weak_candidates']} "
            f"full={a['full_candidates']}"
        )

        shown += 1

        if shown >= 30:
            break

    print()

    # ---------------------------------------------------------
    # Close-factor cases
    # ---------------------------------------------------------

    print("CLOSE-FACTOR CASES")
    print("------------------")

    close_cases = sorted(
        analyses,
        key=lambda a: a["q"] - a["p"]
    )

    for a in close_cases[:30]:
        print(
            f"N={a['N']} "
            f"p={a['p']} "
            f"q={a['q']} "
            f"s={a['s']} "
            f"gap={a['q'] - a['p']} "
            f"q/s={a['q'] / a['s']:.6f} "
            f"t={a['t']} "
            f"x2={a['x2']}"
        )

    print()

    # ---------------------------------------------------------
    # Test simple predictors
    # ---------------------------------------------------------

    print("SIMPLE N/s/D PREDICTORS")
    print("-----------------------")

    predictors = {
        "N mod 2": lambda a: a["N"] % 2,
        "N mod 4": lambda a: a["N"] % 4,
        "s mod 4": lambda a: a["s"] % 4,
        "D mod 4": lambda a: a["D"] % 4,
        "N mod 6": lambda a: a["N"] % 6,
        "D mod 6": lambda a: a["D"] % 6,
        "s mod 6": lambda a: a["s"] % 6,
    }

    for name, predictor in predictors.items():
        values = {}

        for a in analyses:
            key = predictor(a)

            if key not in values:
                values[key] = set()

            values[key].add(a["t"])

        deterministic = sum(
            len(v) == 1
            for v in values.values()
        )

        print(
            f"{name:10s}: "
            f"{deterministic}/{len(values)} groups deterministic"
        )

    print()

    # ---------------------------------------------------------
    # Key exact relationships
    # ---------------------------------------------------------

    print("KEY EXACT RELATIONS")
    print("-------------------")
    print("1. x2 = (t*q - 3*s)/4")
    print("2. 0 <= x2 < q")
    print("3. 3*s <= t*q < 4*q + 3*s")
    print("4. t*q == 3*s (mod 4)")
    print("5. For q>s, the full conditions uniquely determine t.")
    print()

    print("INTERPRETATION")
    print("--------------")
    print("1. The complete branch conditions should select exactly one t.")
    print("2. The interesting problem is whether the necessary q-information")
    print("   can be recovered from N, s, D, or the arithmetic-sequence data.")
    print("3. Pay particular attention to the residue groups and the")
    print("   weakly ambiguous cases.")
    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT_NUMBER}")


if __name__ == "__main__":
    run_experiment()