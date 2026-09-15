import math
import random


EXPERIMENT_NUMBER = 26


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
            is_prime[start:start + count * p:p] = b"\x00" * count

    return [i for i in range(2, limit + 1) if is_prime[i]]


def mod_inverse(a, p):
    """Return inverse of a modulo odd prime p."""
    return pow(a % p, p - 2, p)


def second_branch_root(s, r):
    """
    Solve:

        4*x + 3*s = 0 (mod r)

    using the least nonnegative residue.
    """
    return (-3 * s * mod_inverse(4, r)) % r


def branch_multiplier(s, r, x):
    """
    Given a second-branch root x, compute t from

        4*x + 3*s = t*r.
    """
    numerator = 4 * x + 3 * s

    if numerator % r != 0:
        raise ValueError(
            f"Invalid second-branch root: "
            f"4*{x}+3*{s} is not divisible by {r}"
        )

    return numerator // r


def generate_semiprimes(primes, count, seed):
    """Generate distinct semiprime factor pairs."""
    random.seed(seed)

    pairs = set()

    while len(pairs) < count:
        p = random.choice(primes)
        q = random.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        pairs.add((p, q))

    return sorted(pairs)


def analyze_case(p, q):
    """Analyze second-branch roots for both factors."""
    n = p * q
    s = math.isqrt(n)

    xp2 = second_branch_root(s, p)
    xq2 = second_branch_root(s, q)

    tp = branch_multiplier(s, p, xp2)
    tq = branch_multiplier(s, q, xq2)

    xp = s + 1 - p
    xq = s + 1 - q

    y = p + q - 2 * s - 1

    p_gap = s - p
    q_gap = q - s

    return {
        "n": n,
        "p": p,
        "q": q,
        "s": s,
        "xp2": xp2,
        "xq2": xq2,
        "tp": tp,
        "tq": tq,
        "xp": xp,
        "xq": xq,
        "y": y,
        "p_gap": p_gap,
        "q_gap": q_gap,
        "p_ratio": p / s,
        "q_ratio": q / s,
    }


def histogram(values):
    """Return a frequency dictionary."""
    result = {}

    for value in values:
        result[value] = result.get(value, 0) + 1

    return result


def run_experiment():
    print(f"START EXPERIMENT {EXPERIMENT_NUMBER}")
    print()

    print("SECOND-BRANCH BOUNDED-QUOTIENT EXPERIMENT")
    print("-----------------------------------------")
    print()
    print("For each factor r | N:")
    print("  4*x2 + 3*s = t*r")
    print("  x2 = second-branch root")
    print()
    print("Important correction:")
    print("  The smaller factor p can satisfy p < s.")
    print("  Therefore t is NOT bounded by 6.")
    print()

    primes = sieve_primes(5000)

    # Mix smaller and larger factors.
    cases = generate_semiprimes(
        primes[10:700],
        500,
        seed=26
    )

    analyses = [
        analyze_case(p, q)
        for p, q in cases
    ]

    total = len(analyses)

    tp_values = [a["tp"] for a in analyses]
    tq_values = [a["tq"] for a in analyses]

    tp_hist = histogram(tp_values)
    tq_hist = histogram(tq_values)

    print("MULTIPLIER RANGE")
    print("----------------")

    print(
        f"p branch: min t={min(tp_values)}, "
        f"max t={max(tp_values)}, "
        f"mean t={sum(tp_values) / total:.4f}"
    )

    print(
        f"q branch: min t={min(tq_values)}, "
        f"max t={max(tq_values)}, "
        f"mean t={sum(tq_values) / total:.4f}"
    )

    print()

    print("p-BRANCH t DISTRIBUTION")
    print("-----------------------")

    for t in sorted(tp_hist):
        print(
            f"t={t:4d}: "
            f"{tp_hist[t]:4d}/{total}"
        )

    print()

    print("q-BRANCH t DISTRIBUTION")
    print("-----------------------")

    for t in sorted(tq_hist):
        print(
            f"t={t:4d}: "
            f"{tq_hist[t]:4d}/{total}"
        )

    print()

    print("SECOND-BRANCH ROOT RANGES")
    print("--------------------------")

    p_roots = [a["xp2"] for a in analyses]
    q_roots = [a["xq2"] for a in analyses]

    print(
        f"p branch x2: min={min(p_roots)}, "
        f"max={max(p_roots)}"
    )

    print(
        f"q branch x2: min={min(q_roots)}, "
        f"max={max(q_roots)}"
    )

    print()

    print("t=3 CHECK")
    print("----------")

    p_t3 = [a for a in analyses if a["tp"] == 3]
    q_t3 = [a for a in analyses if a["tq"] == 3]

    print(f"p cases with t=3: {len(p_t3)}")
    print(f"q cases with t=3: {len(q_t3)}")

    p_t3_exact = 0
    for a in p_t3:
        if 4 * a["xp2"] == 3 * (a["p"] - a["s"]):
            p_t3_exact += 1

    q_t3_exact = 0
    for a in q_t3:
        if 4 * a["xq2"] == 3 * (a["q"] - a["s"]):
            q_t3_exact += 1

    print(
        f"p: 4*x2 = 3*(p-s): "
        f"{p_t3_exact}/{len(p_t3)}"
    )

    print(
        f"q: 4*x2 = 3*(q-s): "
        f"{q_t3_exact}/{len(q_t3)}"
    )

    print()

    print("RELATION TO y AND x_p")
    print("---------------------")

    # Exact identity:
    #
    # y = p+q-2s-1
    #
    # x_p = s+1-p
    #
    # therefore:
    #
    # y + x_p = q-s.
    #
    # Thus for t_q=3:
    #
    # 4*x2_q = 3*(y+x_p)

    relation_total = 0
    relation_matches = 0

    for a in analyses:
        if a["tq"] != 3:
            continue

        relation_total += 1

        lhs = 4 * a["xq2"]
        rhs = 3 * (a["y"] + a["xp"])

        if lhs == rhs:
            relation_matches += 1

    print(
        "t_q=3: 4*x2_q = 3*(y+x_p): "
        f"{relation_matches}/{relation_total}"
    )

    print()

    print("t VERSUS FACTOR DISTANCE")
    print("------------------------")

    ranges = [
        ("very close", 0, 50),
        ("close", 50, 200),
        ("medium", 200, 1000),
        ("far", 1000, 10**18),
    ]

    for name, low, high in ranges:
        subset = [
            a
            for a in analyses
            if low <= a["q"] - a["p"] < high
        ]

        if not subset:
            continue

        tq_subset = [a["tq"] for a in subset]

        print(
            f"{name:10s} "
            f"cases={len(subset):4d} "
            f"min_t={min(tq_subset):4d} "
            f"max_t={max(tq_subset):4d} "
            f"mean_t={sum(tq_subset) / len(tq_subset):.4f}"
        )

    print()

    print("CLOSE-FACTOR CASES")
    print("------------------")

    close_cases = sorted(
        analyses,
        key=lambda a: a["q"] - a["p"]
    )

    for a in close_cases[:30]:
        print(
            f"N={a['n']} "
            f"p={a['p']} "
            f"q={a['q']} "
            f"s={a['s']} "
            f"x2p={a['xp2']} "
            f"tp={a['tp']} "
            f"x2q={a['xq2']} "
            f"tq={a['tq']} "
            f"p_gap={a['p_gap']} "
            f"q_gap={a['q_gap']} "
            f"xp={a['xp']} "
            f"y={a['y']}"
        )

    print()

    print("EXTREME MULTIPLIERS")
    print("-------------------")

    extreme = sorted(
        analyses,
        key=lambda a: a["tp"],
        reverse=True
    )

    for a in extreme[:15]:
        print(
            f"N={a['n']} "
            f"p={a['p']} "
            f"q={a['q']} "
            f"s={a['s']} "
            f"x2p={a['xp2']} "
            f"tp={a['tp']}"
        )

    print()

    print("EXACT ALGEBRAIC CHECK")
    print("---------------------")

    # For every case:
    #
    # 4*x2_r + 3*s = t_r*r
    #
    # Verify both branches exactly.
    p_ok = 0
    q_ok = 0

    for a in analyses:
        if (
            4 * a["xp2"] + 3 * a["s"]
            == a["tp"] * a["p"]
        ):
            p_ok += 1

        if (
            4 * a["xq2"] + 3 * a["s"]
            == a["tq"] * a["q"]
        ):
            q_ok += 1

    print(f"p identity checks: {p_ok}/{total}")
    print(f"q identity checks: {q_ok}/{total}")

    print()

    print("INTERPRETATION")
    print("--------------")
    print("1. t is a quotient generated by the second modular branch.")
    print("2. The previous t<=6 assumption was invalid for p<s.")
    print("3. Check whether t_q has a strong pattern for close factors.")
    print("4. Check whether t_q=3 connects x2_q to y+x_p.")
    print("5. The long-term target remains finding x2 or t from N alone.")

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT_NUMBER}")


if __name__ == "__main__":
    run_experiment()