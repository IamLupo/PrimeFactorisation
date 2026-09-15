import math
import random


EXPERIMENT_NUMBER = 25


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
    """Return a^(-1) mod p for prime p."""
    return pow(a, p - 2, p)


def second_branch_residue(s, r):
    """
    Solve:

        4x + 3s = 0 (mod r)

    for x in [0, r-1].
    """
    inv4 = mod_inverse(4 % r, r)
    return (-3 * s * inv4) % r


def first_branch_residue(s, r):
    """
    Solve:

        x = s + 1 (mod r)
    """
    return (s + 1) % r


def factor_coordinate(s, r):
    """
    The coordinate corresponding to factor r:

        x_r = s + 1 - r
    """
    return s + 1 - r


def centered_residue(residue, modulus):
    """
    Return the representative closest to zero.
    """
    if residue > modulus // 2:
        return residue - modulus
    return residue


def binomial_weight(x, y):
    """
    W(x,y) = C(x,y) * 2^(x-y).

    Only defined here on 0 <= y <= x.
    """
    if x < 0 or y < 0 or y > x:
        return None

    return math.comb(x, y) * (1 << (x - y))


def construct_case(p, q):
    """
    Construct all basic quantities.
    """
    n = p * q
    s = math.isqrt(n)
    d = n - s * s
    S = p + q

    x_p = factor_coordinate(s, p)
    x_q = factor_coordinate(s, q)

    y = S - 2 * s - 1

    return {
        "N": n,
        "p": p,
        "q": q,
        "s": s,
        "D": d,
        "S": S,
        "x_p": x_p,
        "x_q": x_q,
        "y": y,
    }


def analyze_branch(case, r):
    """
    Analyze the two modular branches for factor r.
    """
    s = case["s"]

    b1 = first_branch_residue(s, r)
    b2 = second_branch_residue(s, r)

    x_factor = factor_coordinate(s, r)

    # Smallest nonnegative representatives in the full residue class.
    b1_small = b1
    b2_small = b2

    # Distance of the factor coordinate from the second branch.
    # The residue difference is measured modulo r.
    second_delta = (x_factor - b2_small) % r

    # Is the second branch itself inside the natural interval 0..s+1?
    second_in_range = 0 <= b2_small <= s + 1

    # Number of steps of modulus required to reach the actual factor coordinate.
    if (x_factor - b2_small) % r == 0:
        branch_offset = (x_factor - b2_small) // r
    else:
        branch_offset = None

    return {
        "r": r,
        "first": b1_small,
        "second": b2_small,
        "factor_x": x_factor,
        "second_delta_mod_r": second_delta,
        "second_in_range": second_in_range,
        "branch_offset": branch_offset,
    }


def generate_semiprimes(primes, count, seed):
    """
    Generate distinct semiprime factor pairs.
    """
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

    print("H(x) LINEAR-BRANCH EXPERIMENT")
    print("-----------------------------")
    print()
    print("Using:")
    print("  H(x) - 3N = (x-s-1)(4x+3s)")
    print()
    print("For each factor r | N:")
    print("  Branch 1: x = s+1 (mod r)")
    print("  Branch 2: 4x+3s = 0 (mod r)")
    print()
    print("Testing whether Branch 2 has structure connected to:")
    print("  p, q, S=p+q, s, x_p, x_q, and y=S-2s-1")
    print()

    primes = sieve_primes(5000)

    # Avoid too many tiny factors dominating the statistics.
    cases = generate_semiprimes(primes[10:700], 400, seed=25)

    analyses = []

    for p, q in cases:
        case = construct_case(p, q)

        p_branch = analyze_branch(case, p)
        q_branch = analyze_branch(case, q)

        analyses.append(
            {
                "case": case,
                "p_branch": p_branch,
                "q_branch": q_branch,
            }
        )

    total = len(analyses)

    # Aggregate statistics.
    p_second_in_range = sum(
        a["p_branch"]["second_in_range"]
        for a in analyses
    )

    q_second_in_range = sum(
        a["q_branch"]["second_in_range"]
        for a in analyses
    )

    p_second_equals_factor = sum(
        a["p_branch"]["second"] == a["p_branch"]["factor_x"]
        for a in analyses
    )

    q_second_equals_factor = sum(
        a["q_branch"]["second"] == a["q_branch"]["factor_x"]
        for a in analyses
    )

    p_branch_match_xp = sum(
        a["p_branch"]["first"] == a["case"]["x_p"] % a["case"]["p"]
        for a in analyses
    )

    q_branch_match_xq = sum(
        a["q_branch"]["first"] == a["case"]["x_q"] % a["case"]["q"]
        for a in analyses
    )

    print("AGGREGATE RESULTS")
    print("-----------------")
    print(f"Cases tested                         : {total}")
    print(f"p-branch-2 inside 0..s+1            : {p_second_in_range}/{total}")
    print(f"q-branch-2 inside 0..s+1            : {q_second_in_range}/{total}")
    print(f"p branch-2 = actual x_p             : {p_second_equals_factor}/{total}")
    print(f"q branch-2 = actual x_q             : {q_second_equals_factor}/{total}")
    print(f"p first branch matches x_p mod p    : {p_branch_match_xp}/{total}")
    print(f"q first branch matches x_q mod q    : {q_branch_match_xq}/{total}")
    print()

    print("REPRESENTATIVE CASES")
    print("--------------------")

    shown = 0

    for a in analyses:
        c = a["case"]
        pb = a["p_branch"]
        qb = a["q_branch"]

        # Prefer cases where branch 2 is in range for q,
        # since those are the most interesting.
        if not qb["second_in_range"]:
            continue

        print(
            f"N={c['N']} "
            f"p={c['p']} "
            f"q={c['q']} "
            f"s={c['s']} "
            f"S={c['S']} "
            f"y={c['y']}"
        )

        print(
            f"  p: first={pb['first']} "
            f"second={pb['second']} "
            f"x_p={pb['factor_x']} "
            f"offset={pb['branch_offset']}"
        )

        print(
            f"  q: first={qb['first']} "
            f"second={qb['second']} "
            f"x_q={qb['factor_x']} "
            f"offset={qb['branch_offset']}"
        )

        print(
            f"  q second centered="
            f"{centered_residue(qb['second'], c['q'])}"
        )

        print()

        shown += 1

        if shown >= 20:
            break

    print("SECOND-BRANCH RESIDUE STATISTICS")
    print("--------------------------------")
    print()

    # Examine normalized / relative location of the second branch.
    ratios = []

    for a in analyses:
        c = a["case"]
        qb = a["q_branch"]

        if qb["second_in_range"] and c["s"] > 0:
            ratios.append(
                qb["second"] / c["s"]
            )

    if ratios:
        print(f"q second-branch samples : {len(ratios)}")
        print(f"minimum x/s             : {min(ratios):.6f}")
        print(f"maximum x/s             : {max(ratios):.6f}")
        print(f"mean x/s                : {sum(ratios) / len(ratios):.6f}")

    print()

    print("COMPARISON WITH S=p+q")
    print("----------------------")

    # Check simple candidate expressions involving S.
    candidate_names = [
        "S",
        "S-s",
        "S-2s-1",
        "S-s-1",
        "2s+1-S",
    ]

    candidate_match_counts = {
        name: 0
        for name in candidate_names
    }

    q_second_values = []

    for a in analyses:
        c = a["case"]
        x2 = a["q_branch"]["second"]

        if not a["q_branch"]["second_in_range"]:
            continue

        q_second_values.append((c, x2))

        candidates = {
            "S": c["S"],
            "S-s": c["S"] - c["s"],
            "S-2s-1": c["S"] - 2 * c["s"] - 1,
            "S-s-1": c["S"] - c["s"] - 1,
            "2s+1-S": 2 * c["s"] + 1 - c["S"],
        }

        for name, value in candidates.items():
            if x2 == value:
                candidate_match_counts[name] += 1

    for name in candidate_names:
        print(
            f"{name:12s}: "
            f"{candidate_match_counts[name]}/{len(q_second_values)}"
        )

    print()

    print("BINOMIAL / ARITHSEQ CHECK")
    print("-------------------------")
    print()

    binomial_matches = 0
    binomial_total = 0

    for a in analyses:
        c = a["case"]
        x = a["q_branch"]["second"]
        y = c["y"]

        if 0 <= y <= x:
            weight = binomial_weight(x, y)

            if weight is not None:
                # Test whether the second-branch x has
                # a clean divisibility relation with the weight.
                #
                # We test p, q, and N separately.
                checks = (
                    weight % c["p"] == 0,
                    weight % c["q"] == 0,
                    weight % c["N"] == 0,
                )

                binomial_total += 1

                if any(checks):
                    binomial_matches += 1

    print(
        "Second-branch W divisible by p/q/N: "
        f"{binomial_matches}/{binomial_total}"
    )

    print()

    print("MOST INTERESTING q-ROOTS")
    print("-------------------------")

    # Sort by how close the second branch gets to the true factor coordinate.
    interesting = []

    for a in analyses:
        c = a["case"]
        qb = a["q_branch"]

        if qb["second_in_range"]:
            distance = abs(qb["second"] - c["x_p"])

            interesting.append(
                (
                    distance,
                    c,
                    qb,
                )
            )

    interesting.sort(key=lambda item: item[0])

    for distance, c, qb in interesting[:25]:
        print(
            f"N={c['N']} "
            f"p={c['p']} "
            f"q={c['q']} "
            f"s={c['s']} "
            f"x_p={c['x_p']} "
            f"x2_q={qb['second']} "
            f"|difference|={distance} "
            f"y={c['y']}"
        )

    print()
    print("KEY EXACT IDENTITIES")
    print("--------------------")
    print(
        "1. H(x)-3N = (x-s-1)(4x+3s)"
    )
    print(
        "2. x_p = s+1-p"
    )
    print(
        "3. x_p is always in Branch 1 modulo p"
    )
    print(
        "4. Branch 2 modulo r is x = -3s * 4^(-1) (mod r)"
    )
    print(
        "5. The experiment asks whether Branch 2 carries"
        " information about the hidden factors beyond the trivial congruence."
    )

    print()
    print("FINISHED EXPERIMENT 25")


if __name__ == "__main__":
    run_experiment()
