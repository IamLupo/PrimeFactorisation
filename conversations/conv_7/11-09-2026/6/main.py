import math
import random


EXPERIMENT_NUMBER = 28


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


def second_branch_root(s, r):
    """Solve 4*x + 3*s = 0 (mod r)."""
    return (-3 * s * pow(4, r - 2, r)) % r


def branch_multiplier(s, r, x):
    """Return t from 4*x + 3*s = t*r."""
    value = 4 * x + 3 * s

    if value % r != 0:
        raise ValueError("Invalid second-branch root.")

    return value // r


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
    """Construct the full coordinate system."""
    n = p * q
    s = math.isqrt(n)

    xp = s + 1 - p
    y = p + q - 2 * s - 1

    x2 = second_branch_root(s, q)
    t = branch_multiplier(s, q, x2)

    k = s + 1 - x2

    lhs = 4 * (k - p)

    rhs = (3 - t) * s + (4 - t) * xp - t * y

    return {
        "n": n,
        "p": p,
        "q": q,
        "s": s,
        "xp": xp,
        "y": y,
        "x2": x2,
        "t": t,
        "k": k,
        "k_minus_p": k - p,
        "lhs": lhs,
        "rhs": rhs,
    }


def run_experiment():
    print(f"START EXPERIMENT {EXPERIMENT_NUMBER}")
    print()

    print("SECOND-BRANCH COORDINATE IDENTITY EXPERIMENT")
    print("---------------------------------------------")
    print()
    print("Definitions:")
    print("  x_p = s + 1 - p")
    print("  y   = p + q - 2s - 1")
    print("  4*x2 + 3s = t*q")
    print("  k   = s + 1 - x2")
    print()
    print("Target identity:")
    print("  4(k-p) = (3-t)s + (4-t)x_p - t*y")
    print()
    print("Special t=3 identity:")
    print("  4(k-p) = x_p - 3y")
    print()

    primes = sieve_primes(5000)

    cases = generate_semiprimes(
        primes[10:700],
        1000,
        seed=28
    )

    analyses = [
        analyze_case(p, q)
        for p, q in cases
    ]

    total = len(analyses)

    print("AGGREGATE RESULTS")
    print("-----------------")

    identity_matches = sum(
        a["lhs"] == a["rhs"]
        for a in analyses
    )

    print(
        f"General identity verified        : "
        f"{identity_matches}/{total}"
    )

    print()

    print("BY t")
    print("----")

    for t in sorted({a["t"] for a in analyses}):
        subset = [a for a in analyses if a["t"] == t]

        matches = sum(
            a["lhs"] == a["rhs"]
            for a in subset
        )

        print(
            f"t={t}: "
            f"cases={len(subset)} "
            f"identity={matches}/{len(subset)}"
        )

    print()

    print("t=3 SPECIAL IDENTITY")
    print("--------------------")

    t3 = [a for a in analyses if a["t"] == 3]

    t3_matches = 0

    for a in t3:
        if 4 * (a["k"] - a["p"]) == a["xp"] - 3 * a["y"]:
            t3_matches += 1

    print(f"t=3 cases                   : {len(t3)}")
    print(
        f"4(k-p)=x_p-3y               : "
        f"{t3_matches}/{len(t3)}"
    )

    print()

    print("H(x_p)=0 CONDITION")
    print("-------------------")

    zero_h = [
        a for a in analyses
        if a["k"] == a["p"]
    ]

    print(
        f"k=p cases                   : "
        f"{len(zero_h)}/{total}"
    )

    # For t=3:
    # k=p <=> x_p=3y.
    t3_k_equals_p = [
        a for a in t3
        if a["k"] == a["p"]
    ]

    t3_xp_equals_3y = [
        a for a in t3
        if a["xp"] == 3 * a["y"]
    ]

    print(
        f"t=3 and k=p                 : "
        f"{len(t3_k_equals_p)}/{len(t3)}"
    )

    print(
        f"t=3 and x_p=3y              : "
        f"{len(t3_xp_equals_3y)}/{len(t3)}"
    )

    print()

    print("ALL t=3 CASES WITH x_p=3y")
    print("-------------------------")

    for a in t3_xp_equals_3y:
        print(
            f"N={a['n']} "
            f"p={a['p']} "
            f"q={a['q']} "
            f"s={a['s']} "
            f"x_p={a['xp']} "
            f"y={a['y']} "
            f"x2={a['x2']} "
            f"k={a['k']} "
            f"t={a['t']}"
        )

    print()

    print("LINEAR FORMS BY t")
    print("-----------------")
    print()
    print("For each t:")
    print("  L_t = (3-t)s + (4-t)x_p - t*y")
    print("and test whether L_t has a simpler factor-distance relation.")

    for t in sorted({a["t"] for a in analyses}):
        subset = [a for a in analyses if a["t"] == t]

        if not subset:
            continue

        print()
        print(f"t={t}")

        # Show the first few exact values.
        for a in subset[:8]:
            print(
                f"  N={a['n']} "
                f"s={a['s']} "
                f"x_p={a['xp']} "
                f"y={a['y']} "
                f"k-p={a['k_minus_p']} "
                f"L={a['rhs']}"
            )

    print()

    print("CLOSE-FACTOR t=3 CASES")
    print("-----------------------")

    close_t3 = sorted(
        t3,
        key=lambda a: a["q"] - a["p"]
    )

    for a in close_t3[:30]:
        print(
            f"N={a['n']} "
            f"p={a['p']} "
            f"q={a['q']} "
            f"s={a['s']} "
            f"x_p={a['xp']} "
            f"y={a['y']} "
            f"x2={a['x2']} "
            f"k-p={a['k_minus_p']}"
        )

    print()

    print("SMALLEST |k-p|")
    print("--------------")

    closest = sorted(
        analyses,
        key=lambda a: abs(a["k_minus_p"])
    )

    for a in closest[:30]:
        print(
            f"N={a['n']} "
            f"p={a['p']} "
            f"q={a['q']} "
            f"s={a['s']} "
            f"t={a['t']} "
            f"x_p={a['xp']} "
            f"y={a['y']} "
            f"k-p={a['k_minus_p']}"
        )

    print()

    print("INTERPRETATION")
    print("--------------")
    print("1. Verify the general linear identity for every t.")
    print("2. Determine whether t=3 gives a special coordinate relation.")
    print("3. Test whether H(x_p)=0 is equivalent to x_p=3y inside t=3.")
    print("4. Examine whether other t values yield similarly simple relations.")
    print("5. Look for a relation that determines x_p from s and y.")

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT_NUMBER}")


if __name__ == "__main__":
    run_experiment()
