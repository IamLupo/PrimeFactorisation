import math
import random


EXPERIMENT_NUMBER = 27


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
            is_prime[start:start + count * p:start if False else p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i in range(2, limit + 1) if is_prime[i]]


def second_branch_root(s, r):
    """
    Solve:

        4*x + 3*s = 0 (mod r)

    using the least nonnegative residue.
    """
    return (-3 * s * pow(4, r - 2, r)) % r


def branch_multiplier(s, r, x):
    """
    Compute t from:

        4*x + 3*s = t*r
    """
    value = 4 * x + 3 * s

    if value % r != 0:
        raise ValueError("Invalid second-branch root.")

    return value // r


def construct_case(p, q):
    """Construct all factor coordinates."""
    n = p * q
    s = math.isqrt(n)

    x_p = s + 1 - p
    y = p + q - 2 * s - 1

    x2_q = second_branch_root(s, q)
    t_q = branch_multiplier(s, q, x2_q)

    k_q = s + 1 - x2_q

    # H evaluated at the true p-coordinate.
    H_xp = (
        4 * x_p * x_p
        - (s + 4) * x_p
        + 3 * (n - s * s)
        - 3 * s
    )

    return {
        "N": n,
        "p": p,
        "q": q,
        "s": s,
        "x_p": x_p,
        "y": y,
        "x2_q": x2_q,
        "t_q": t_q,
        "k_q": k_q,
        "H_xp": H_xp,
    }


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


def run_experiment():
    print(f"START EXPERIMENT {EXPERIMENT_NUMBER}")
    print()

    print("SECOND-BRANCH / FACTOR-COORDINATE EXPERIMENT")
    print("---------------------------------------------")
    print()
    print("Definitions:")
    print("  x_p  = s + 1 - p")
    print("  4*x2 + 3*s = t*q")
    print("  k    = s + 1 - x2")
    print()
    print("Then:")
    print("  4*(k-p) = 7s + 4 - t*q - 4p")
    print()
    print("For t=3:")
    print("  H(x_p) = -4*p*(k-p)")
    print()
    print("The experiment tests whether k-p, t, y, and H(x_p)")
    print("show a systematic relationship.")
    print()

    primes = sieve_primes(5000)

    cases = generate_semiprimes(
        primes[10:700],
        500,
        seed=27
    )

    analyses = [
        construct_case(p, q)
        for p, q in cases
    ]

    total = len(analyses)

    print("AGGREGATE RESULTS")
    print("-----------------")

    t3_cases = [
        a for a in analyses
        if a["t_q"] == 3
    ]

    print(f"Cases tested                  : {total}")
    print(f"q-branch t=3                 : {len(t3_cases)}/{total}")

    # Exact t=3 relation.
    relation_ok = 0

    for a in t3_cases:
        lhs = 4 * a["x2_q"]
        rhs = 3 * (a["q"] - a["s"])

        if lhs == rhs:
            relation_ok += 1

    print(
        f"4*x2 = 3*(q-s)               : "
        f"{relation_ok}/{len(t3_cases)}"
    )

    # Main new identity.
    identity_ok = 0

    for a in t3_cases:
        p = a["p"]
        k = a["k_q"]

        lhs = a["H_xp"]
        rhs = -4 * p * (k - p)

        if lhs == rhs:
            identity_ok += 1

    print(
        f"H(x_p)=-4p(k-p)              : "
        f"{identity_ok}/{len(t3_cases)}"
    )

    # Check y + xp = q-s.
    y_relation_ok = 0

    for a in analyses:
        lhs = a["y"] + a["x_p"]
        rhs = a["q"] - a["s"]

        if lhs == rhs:
            y_relation_ok += 1

    print(
        f"y+x_p=q-s                    : "
        f"{y_relation_ok}/{total}"
    )

    print()

    print("K-P DISTRIBUTION FOR t=3")
    print("------------------------")

    differences = {}

    for a in t3_cases:
        difference = a["k_q"] - a["p"]
        differences[difference] = differences.get(difference, 0) + 1

    for difference in sorted(differences):
        print(
            f"k-p={difference:6d}: "
            f"{differences[difference]:4d}/{len(t3_cases)}"
        )

    print()

    print("RELATION TO y")
    print("-------------")

    # For t=3:
    #
    # x2 = 3(q-s)/4
    #
    # q-s = y+x_p
    #
    # k = s+1-x2
    #
    # Therefore:
    #
    # k = s+1 - 3(y+x_p)/4.
    #
    # Test this exact expression.
    y_formula_ok = 0

    for a in t3_cases:
        if (
            4 * a["k_q"]
            == 4 * (a["s"] + 1)
            - 3 * (a["y"] + a["x_p"])
        ):
            y_formula_ok += 1

    print(
        f"4k = 4(s+1)-3(y+x_p)        : "
        f"{y_formula_ok}/{len(t3_cases)}"
    )

    print()

    print("CLOSE-FACTOR CASES")
    print("------------------")

    close_cases = sorted(
        analyses,
        key=lambda a: a["q"] - a["p"]
    )

    for a in close_cases[:40]:
        k_minus_p = a["k_q"] - a["p"]

        print(
            f"N={a['N']} "
            f"p={a['p']} "
            f"q={a['q']} "
            f"s={a['s']} "
            f"t={a['t_q']} "
            f"x_p={a['x_p']} "
            f"x2={a['x2_q']} "
            f"k={a['k_q']} "
            f"k-p={k_minus_p} "
            f"y={a['y']} "
            f"Hxp={a['H_xp']}"
        )

    print()

    print("CASES WHERE k = p")
    print("-----------------")

    k_equals_p = [
        a for a in analyses
        if a["k_q"] == a["p"]
    ]

    print(
        f"k=p cases: {len(k_equals_p)}/{total}"
    )

    for a in k_equals_p[:30]:
        print(
            f"N={a['N']} "
            f"p={a['p']} "
            f"q={a['q']} "
            f"s={a['s']} "
            f"t={a['t_q']} "
            f"x_p={a['x_p']} "
            f"x2={a['x2_q']} "
            f"y={a['y']} "
            f"Hxp={a['H_xp']}"
        )

    print()

    print("CASES WHERE H(x_p)=0")
    print("--------------------")

    h_zero = [
        a for a in analyses
        if a["H_xp"] == 0
    ]

    print(
        f"H(x_p)=0 cases: {len(h_zero)}/{total}"
    )

    for a in h_zero[:30]:
        print(
            f"N={a['N']} "
            f"p={a['p']} "
            f"q={a['q']} "
            f"s={a['s']} "
            f"t={a['t_q']} "
            f"x_p={a['x_p']} "
            f"x2={a['x2_q']} "
            f"k={a['k_q']} "
            f"k-p={a['k_q'] - a['p']} "
            f"y={a['y']}"
        )

    print()

    print("SECOND-BRANCH CONSTRUCTION CHECK")
    print("--------------------------------")

    construction_ok = 0

    for a in analyses:
        if (
            4 * a["x2_q"] + 3 * a["s"]
            == a["t_q"] * a["q"]
        ):
            construction_ok += 1

    print(
        f"4*x2+3s=tq: "
        f"{construction_ok}/{total}"
    )

    print()

    print("INTERPRETATION")
    print("--------------")
    print("1. Determine whether k=s+1-x2 has structure closer to p.")
    print("2. For t=3, test the exact H(x_p)=-4p(k-p) identity.")
    print("3. Check whether k=p is associated with H(x_p)=0.")
    print("4. Examine whether k-p is related to y or factor distance.")
    print("5. The important question remains whether k can be predicted from N alone.")

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT_NUMBER}")


if __name__ == "__main__":
    run_experiment()
