import math
import random


EXPERIMENT_NUMBER = 31


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
    """Solve 4*x + 3*s = 0 (mod q)."""
    return (-3 * s * pow(4, q - 2, q)) % q


def branch_multiplier(s, q, x):
    """Compute t from 4*x + 3*s = t*q."""
    value = 4 * x + 3 * s

    if value % q != 0:
        raise ValueError("Invalid second-branch root.")

    return value // q


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
    """Construct the quantities used by the experiment."""
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
        "q_mod_4": q % 4,
    }


def deterministic_group_statistics(analyses, modulus):
    """
    Group cases by N mod M, s mod M, D mod M.

    Return:
        number of groups
        number of deterministic groups
        number of cases belonging to deterministic groups
        number of total cases
    """
    groups = {}

    for a in analyses:
        key = (
            a["N"] % modulus,
            a["s"] % modulus,
            a["D"] % modulus,
        )

        groups.setdefault(key, set()).add(a["t"])

    deterministic_groups = 0
    deterministic_cases = 0

    for key, values in groups.items():
        if len(values) == 1:
            deterministic_groups += 1

            # We need the number of cases in the group,
            # so this is computed separately below.
            #
            # Kept here only for the group count.
            continue

    case_counts = {}

    for a in analyses:
        key = (
            a["N"] % modulus,
            a["s"] % modulus,
            a["D"] % modulus,
        )

        case_counts[key] = case_counts.get(key, 0) + 1

    for key, values in groups.items():
        if len(values) == 1:
            deterministic_cases += case_counts[key]

    return (
        len(groups),
        deterministic_groups,
        deterministic_cases,
        len(analyses),
    )


def group_t_distribution(analyses, modulus):
    """
    Return groups whose same (N,s,D) residue fingerprint
    contains multiple t values.
    """
    groups = {}

    for a in analyses:
        key = (
            a["N"] % modulus,
            a["s"] % modulus,
            a["D"] % modulus,
        )

        groups.setdefault(key, set()).add(a["t"])

    return groups


def run_experiment():
    print(f"START EXPERIMENT {EXPERIMENT_NUMBER}")
    print()

    print("2-ADIC BRANCH-FINGERPRINT EXPERIMENT")
    print("------------------------------------")
    print()
    print("For each semiprime:")
    print("  N = p*q")
    print("  s = floor(sqrt(N))")
    print("  D = N-s^2")
    print("  4*x2+3*s = t*q")
    print()
    print("Question:")
    print("  Does increasing 2-adic information from")
    print("  (N,s,D) make t increasingly deterministic?")
    print()

    primes = sieve_primes(5000)

    cases = generate_semiprimes(
        primes,
        2000,
        seed=31
    )

    analyses = [
        analyze_case(p, q)
        for p, q in cases
    ]

    total = len(analyses)

    print("AGGREGATE RESULTS")
    print("-----------------")
    print(f"Cases tested: {total}")
    print()

    # ---------------------------------------------------------
    # Test powers of two.
    # ---------------------------------------------------------

    print("2-ADIC FINGERPRINT")
    print("------------------")

    moduli = [
        2,
        4,
        8,
        16,
        32,
        64,
        128,
        256,
        512,
        1024,
        2048,
        4096,
    ]

    for modulus in moduli:
        (
            group_count,
            deterministic_group_count,
            deterministic_case_count,
            case_count,
        ) = deterministic_group_statistics(
            analyses,
            modulus
        )

        group_percentage = (
            100.0
            * deterministic_group_count
            / group_count
        )

        case_percentage = (
            100.0
            * deterministic_case_count
            / case_count
        )

        print(
            f"M={modulus:4d} "
            f"groups={group_count:4d} "
            f"det_groups={deterministic_group_count:4d} "
            f"det_group%={group_percentage:7.3f} "
            f"det_cases={deterministic_case_count:4d}/{case_count} "
            f"det_case%={case_percentage:7.3f}"
        )

    print()

    # ---------------------------------------------------------
    # Compare two fingerprints:
    #
    # A = (N,s)
    # B = (N,s,D)
    #
    # D is algebraically determined by N and s, but keeping
    # it explicit makes the result easier to inspect.
    # ---------------------------------------------------------

    print("COMPARISON: (N,s) VS (N,s,D)")
    print("---------------------------")

    for modulus in [4, 8, 16, 32, 64, 128, 256, 512]:
        groups_ns = {}
        groups_nsd = {}

        for a in analyses:
            key_ns = (
                a["N"] % modulus,
                a["s"] % modulus,
            )

            key_nsd = (
                a["N"] % modulus,
                a["s"] % modulus,
                a["D"] % modulus,
            )

            groups_ns.setdefault(key_ns, set()).add(a["t"])
            groups_nsd.setdefault(key_nsd, set()).add(a["t"])

        ns_det = sum(
            len(values) == 1
            for values in groups_ns.values()
        )

        nsd_det = sum(
            len(values) == 1
            for values in groups_nsd.values()
        )

        print(
            f"M={modulus:4d} "
            f"(N,s) det={ns_det:4d}/{len(groups_ns):4d} "
            f"(N,s,D) det={nsd_det:4d}/{len(groups_nsd):4d}"
        )

    print()

    # ---------------------------------------------------------
    # Find the first modulus where every observed group
    # becomes deterministic.
    # ---------------------------------------------------------

    print("FIRST FULLY-DETERMINISTIC MODULUS")
    print("---------------------------------")

    first_full = None

    for exponent in range(1, 21):
        modulus = 2 ** exponent

        groups = group_t_distribution(
            analyses,
            modulus
        )

        if all(len(values) == 1 for values in groups.values()):
            first_full = modulus
            break

    if first_full is None:
        print(
            "No modulus 2^k with k<=20 made every observed "
            "fingerprint deterministic."
        )
    else:
        print(
            f"First fully deterministic modulus: {first_full}"
        )

    print()

    # ---------------------------------------------------------
    # Show the actual ambiguous groups at M=16 and M=256.
    # ---------------------------------------------------------

    for modulus in [16, 256]:
        print(
            f"AMBIGUOUS GROUPS AT M={modulus}"
        )
        print(
            "-" * (29 + len(str(modulus)))
        )

        groups = group_t_distribution(
            analyses,
            modulus
        )

        ambiguous = [
            (key, values)
            for key, values in groups.items()
            if len(values) > 1
        ]

        ambiguous.sort(
            key=lambda item: (
                item[0],
                sorted(item[1])
            )
        )

        for key, values in ambiguous[:30]:
            print(
                f"fingerprint={key} "
                f"possible_t={sorted(values)}"
            )

        print(
            f"Total ambiguous groups: "
            f"{len(ambiguous)}"
        )

        print()

    # ---------------------------------------------------------
    # q mod 4 versus 2-adic fingerprint.
    # ---------------------------------------------------------

    print("q MOD 4 RECONSTRUCTION TEST")
    print("----------------------------")

    for modulus in [4, 8, 16, 32, 64, 128, 256]:
        groups = {}

        for a in analyses:
            key = (
                a["N"] % modulus,
                a["s"] % modulus,
                a["D"] % modulus,
            )

            groups.setdefault(key, set()).add(
                a["q_mod_4"]
            )

        deterministic = sum(
            len(values) == 1
            for values in groups.values()
        )

        print(
            f"M={modulus:4d} "
            f"qmod4 deterministic groups="
            f"{deterministic}/{len(groups)}"
        )

    print()

    # ---------------------------------------------------------
    # Representative ambiguous examples.
    # ---------------------------------------------------------

    print("REPRESENTATIVE AMBIGUOUS CASES")
    print("------------------------------")

    modulus = 256

    groups = {}

    for index, a in enumerate(analyses):
        key = (
            a["N"] % modulus,
            a["s"] % modulus,
            a["D"] % modulus,
        )

        groups.setdefault(key, []).append(index)

    shown = 0

    for key, indices in groups.items():
        t_values = sorted({
            analyses[i]["t"]
            for i in indices
        })

        if len(t_values) <= 1:
            continue

        print(
            f"fingerprint={key} "
            f"possible_t={t_values}"
        )

        for i in indices[:6]:
            a = analyses[i]

            print(
                f"  N={a['N']} "
                f"p={a['p']} "
                f"q={a['q']} "
                f"s={a['s']} "
                f"t={a['t']} "
                f"qmod4={a['q_mod_4']}"
            )

        print()

        shown += 1

        if shown >= 15:
            break

    # ---------------------------------------------------------
    # t distribution by q mod 4.
    # ---------------------------------------------------------

    print("t DISTRIBUTION BY q MOD 4")
    print("-------------------------")

    for q_mod4 in [1, 3]:
        histogram = {}

        for a in analyses:
            if a["q_mod_4"] != q_mod4:
                continue

            histogram[a["t"]] = (
                histogram.get(a["t"], 0) + 1
            )

        distribution = " ".join(
            f"t{t}={histogram[t]}"
            for t in sorted(histogram)
        )

        count = sum(histogram.values())

        print(
            f"q mod4={q_mod4} "
            f"cases={count} "
            f"{distribution}"
        )

    print()

    print("KEY OBSERVATION TO WATCH")
    print("------------------------")
    print("If larger powers of two substantially increase")
    print("the deterministic-group percentage, that suggests")
    print("the hidden branch is leaving a 2-adic fingerprint.")
    print()
    print("If the percentage quickly plateaus, then the missing")
    print("q mod 4 information is probably not recoverable from")
    print("N,s,D through simple 2-adic data alone.")
    print()

    print("INTERPRETATION")
    print("--------------")
    print("1. t is already exactly characterized if q is known.")
    print("2. This experiment tests whether N and sqrt(N)-derived")
    print("   residues leak enough information to identify t.")
    print("3. Pay attention to the transition M=4 -> 8 -> 16 -> ...")
    print("4. Pay attention to whether q mod 4 becomes deterministic.")
    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT_NUMBER}")


if __name__ == "__main__":
    run_experiment()
