import math
import random


EXPERIMENT_NUMBER = 29


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
            is_prime[start:start + count * p:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i in range(2, limit + 1) if is_prime[i]]


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
        raise ValueError("x is not a valid second-branch root.")

    return value // q


def predicted_t_candidates_from_q(s, q):
    """
    Return all t in {1,...,6} satisfying both:

        t*q >= 3*s
        t*q ≡ 3*s (mod 4)

    This is the complete branch prediction when q is known.
    """
    candidates = []

    for t in range(1, 7):
        if t * q < 3 * s:
            continue

        if (t * q - 3 * s) % 4 != 0:
            continue

        candidates.append(t)

    return candidates


def parity_candidate_class(s, q):
    """
    Return t mod 4 implied by:

        t*q ≡ 3*s (mod 4)

    Since q is odd, q^{-1} ≡ q (mod 4).
    """
    return (3 * (s % 4) * (q % 4)) % 4


def q_ratio_lower_bound(s, q):
    """
    ceil(3*s/q) is the smallest integer satisfying t*q >= 3*s.
    """
    return (3 * s + q - 1) // q


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

    candidates = predicted_t_candidates_from_q(s, q)
    lower = q_ratio_lower_bound(s, q)
    residue_class = parity_candidate_class(s, q)

    return {
        "N": n,
        "p": p,
        "q": q,
        "s": s,
        "D": d,
        "x2": x2,
        "t": t,
        "candidates": candidates,
        "lower": lower,
        "residue_class": residue_class,
        "N_mod4": n % 4,
        "s_mod4": s % 4,
        "D_mod4": d % 4,
        "q_mod4": q % 4,
    }


def run_experiment():
    print(f"START EXPERIMENT {EXPERIMENT_NUMBER}")
    print()

    print("SECOND-BRANCH t CLASSIFICATION EXPERIMENT")
    print("-----------------------------------------")
    print()
    print("Second branch:")
    print("  4*x2 + 3*s = t*q")
    print()
    print("Hypothesis:")
    print("  1 <= t <= 6")
    print("  t*q >= 3*s")
    print("  t*q == 3*s (mod 4)")
    print()
    print("Therefore t should be the unique candidate in {1,...,6}.")
    print()

    primes = sieve_primes(5000)

    cases = generate_semiprimes(
        primes,
        1500,
        seed=29
    )

    analyses = [
        analyze_case(p, q)
        for p, q in cases
    ]

    total = len(analyses)

    # ---------------------------------------------------------
    # Exact candidate prediction
    # ---------------------------------------------------------

    unique_candidate = 0
    actual_in_candidates = 0

    for a in analyses:
        candidates = a["candidates"]

        if len(candidates) == 1:
            unique_candidate += 1

        if a["t"] in candidates:
            actual_in_candidates += 1

    print("AGGREGATE RESULTS")
    print("-----------------")
    print(f"Cases tested                         : {total}")
    print(
        f"Actual t in predicted candidates    : "
        f"{actual_in_candidates}/{total}"
    )
    print(
        f"Exactly one candidate               : "
        f"{unique_candidate}/{total}"
    )
    print()

    # ---------------------------------------------------------
    # Compare actual t to ceil(3s/q)
    # ---------------------------------------------------------

    print("LOWER-BOUND RELATION")
    print("--------------------")

    lower_matches = 0

    for a in analyses:
        if a["t"] == a["lower"] or a["t"] > a["lower"]:
            lower_matches += 1

    print(
        f"t >= ceil(3s/q)                    : "
        f"{lower_matches}/{total}"
    )

    print()

    # ---------------------------------------------------------
    # t distribution
    # ---------------------------------------------------------

    print("ACTUAL t DISTRIBUTION")
    print("---------------------")

    histogram = {}

    for a in analyses:
        histogram[a["t"]] = histogram.get(a["t"], 0) + 1

    for t in sorted(histogram):
        print(
            f"t={t}: {histogram[t]:4d}/{total}"
        )

    print()

    # ---------------------------------------------------------
    # Candidate multiplicity
    # ---------------------------------------------------------

    print("NUMBER OF POSSIBLE t VALUES")
    print("----------------------------")

    candidate_histogram = {}

    for a in analyses:
        count = len(a["candidates"])
        candidate_histogram[count] = (
            candidate_histogram.get(count, 0) + 1
        )

    for count in sorted(candidate_histogram):
        print(
            f"{count} candidate(s): "
            f"{candidate_histogram[count]:4d}/{total}"
        )

    print()

    # ---------------------------------------------------------
    # Explicit residue analysis
    # ---------------------------------------------------------

    print("t MOD 4 ANALYSIS")
    print("Checking:")
    print("  t*q = 3*s (mod 4)")
    print()

    residue_ok = 0

    for a in analyses:
        if (
            a["t"] * a["q"] - 3 * a["s"]
        ) % 4 == 0:
            residue_ok += 1

    print(
        f"Congruence verified                    : "
        f"{residue_ok}/{total}"
    )

    print()

    # ---------------------------------------------------------
    # Try to infer possible t-class from N and s alone.
    #
    # For odd p,q:
    #   N mod 4 = pq mod 4
    #
    # But N mod 4 alone does not generally recover q mod 4.
    # We therefore explicitly test how much information remains.
    # ---------------------------------------------------------

    print("N MOD 4 / S MOD 4 GROUPS")
    print("------------------------")

    groups = {}

    for a in analyses:
        key = (a["N_mod4"], a["s_mod4"])

        if key not in groups:
            groups[key] = []

        groups[key].append(a["t"])

    for key in sorted(groups):
        values = groups[key]

        local_hist = {}
        for t in values:
            local_hist[t] = local_hist.get(t, 0) + 1

        distribution = " ".join(
            f"t{t}={local_hist[t]}"
            for t in sorted(local_hist)
        )

        print(
            f"N mod4={key[0]} "
            f"s mod4={key[1]} "
            f"cases={len(values)} "
            f"{distribution}"
        )

    print()

    # ---------------------------------------------------------
    # D mod 4
    # ---------------------------------------------------------

    print("D MOD 4 GROUPS")
    print("--------------")

    d_groups = {}

    for a in analyses:
        key = (a["s_mod4"], a["D_mod4"])

        if key not in d_groups:
            d_groups[key] = []

        d_groups[key].append(a["t"])

    for key in sorted(d_groups):
        values = d_groups[key]

        local_hist = {}
        for t in values:
            local_hist[t] = local_hist.get(t, 0) + 1

        distribution = " ".join(
            f"t{t}={local_hist[t]}"
            for t in sorted(local_hist)
        )

        print(
            f"s mod4={key[0]} "
            f"D mod4={key[1]} "
            f"cases={len(values)} "
            f"{distribution}"
        )

    print()

    # ---------------------------------------------------------
    # Representative cases
    # ---------------------------------------------------------

    print("REPRESENTATIVE CASES")
    print("--------------------")

    # Show examples covering all t values.
    shown_t = set()

    for a in analyses:
        t = a["t"]

        if t in shown_t:
            continue

        print(
            f"N={a['N']} "
            f"p={a['p']} "
            f"q={a['q']} "
            f"s={a['s']} "
            f"D={a['D']} "
            f"x2={a['x2']} "
            f"t={a['t']} "
            f"lower={a['lower']} "
            f"candidates={a['candidates']} "
            f"Nmod4={a['N_mod4']} "
            f"smod4={a['s_mod4']} "
            f"Dmod4={a['D_mod4']}"
        )

        shown_t.add(t)

        if len(shown_t) == 6:
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
            f"t={a['t']} "
            f"x2={a['x2']} "
            f"lower={a['lower']} "
            f"candidates={a['candidates']}"
        )

    print()

    # ---------------------------------------------------------
    # Search for simple N-only predictors.
    #
    # We deliberately test only small deterministic expressions.
    # ---------------------------------------------------------

    print("SIMPLE N-ONLY CANDIDATE TESTS")
    print("-----------------------------")

    candidate_predictors = {
        "1 + (N mod 6)": lambda a: 1 + (a["N"] % 6),
        "1 + (s mod 6)": lambda a: 1 + (a["s"] % 6),
        "1 + (D mod 6)": lambda a: 1 + (a["D"] % 6),
        "1 + (N mod 4)": lambda a: 1 + (a["N"] % 4),
        "1 + (s mod 4)": lambda a: 1 + (a["s"] % 4),
        "1 + (D mod 4)": lambda a: 1 + (a["D"] % 4),
    }

    for name, predictor in candidate_predictors.items():
        matches = 0
        valid = 0

        for a in analyses:
            predicted = predictor(a)

            if predicted < 1 or predicted > 6:
                continue

            valid += 1

            if predicted == a["t"]:
                matches += 1

        print(
            f"{name:20s}: "
            f"{matches}/{valid} "
            f"valid"
        )

    print()

    print("INTERPRETATION")
    print("--------------")
    print("1. Check whether the range + mod-4 conditions uniquely determine t.")
    print("2. If they do, t is fully characterized once q is known.")
    print("3. Check whether N, s, and D modulo small numbers predict t without q.")
    print("4. A failure of N-only predictors is useful: it tells us the hidden q")
    print("   residue is not trivially encoded by these small residues.")
    print("5. The next useful direction would be replacing q mod 4 with a quantity")
    print("   computable from N and the arithmetic-sequence coordinates.")

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT_NUMBER}")


if __name__ == "__main__":
    run_experiment()
