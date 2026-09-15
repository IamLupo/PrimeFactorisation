import math
import random
from collections import defaultdict


EXPERIMENT_NUMBER = 32


def is_probable_prime(n):
    """
    Deterministic Miller-Rabin for unsigned 64-bit integers.
    """
    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
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

    # Deterministic for n < 2^64.
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
    """
    Generate a random prime with the requested bit length.
    """
    if bits < 3:
        raise ValueError("bits must be >= 3")

    while True:
        candidate = rng.getrandbits(bits)

        # Force exact bit length and oddness.
        candidate |= (1 << (bits - 1))
        candidate |= 1

        if is_probable_prime(candidate):
            return candidate


def generate_semiprime(bits_p, bits_q, rng):
    """
    Generate distinct primes p < q with approximately
    the requested bit sizes.
    """
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

    using the least nonnegative residue.
    """
    return (-3 * s * pow(4, q - 2, q)) % q


def branch_multiplier(s, q, x2):
    """
    Compute t from:

        4*x2 + 3*s = t*q
    """
    value = 4 * x2 + 3 * s

    if value % q != 0:
        raise ValueError(
            "Invalid second-branch root."
        )

    return value // q


def construct_case(p, q):
    """
    Construct the quantities needed for one sample.
    """
    n = p * q
    s = math.isqrt(n)

    x2 = second_branch_root(s, q)
    t = branch_multiplier(s, q, x2)

    return {
        "N": n,
        "p": p,
        "q": q,
        "s": s,
        "t": t,
        "bits_N": n.bit_length(),
    }


def fingerprint(case, k):
    """
    Return the k-bit fingerprint:

        (N mod 2^k, s mod 2^k)
    """
    mask = (1 << k) - 1

    return (
        case["N"] & mask,
        case["s"] & mask,
    )


def fingerprint_is_deterministic(cases, k):
    """
    Determine whether the fingerprint uniquely determines t
    inside the supplied dataset.

    Returns:
        deterministic,
        groups
    """
    groups = {}

    for case in cases:
        key = fingerprint(case, k)
        t = case["t"]

        if key not in groups:
            groups[key] = t
        elif groups[key] != t:
            return False, groups

    return True, groups


def minimal_k_per_case(cases, max_k):
    """
    For each case, find the smallest k such that no other
    case in the dataset shares its k-bit fingerprint with
    a different t.

    Returns:
        dictionary indexed by case index.
    """
    unresolved = set(range(len(cases)))
    result = {}

    for k in range(1, max_k + 1):
        groups = defaultdict(list)

        for index in unresolved:
            key = fingerprint(cases[index], k)
            groups[key].append(index)

        next_unresolved = set()

        for indices in groups.values():
            t_values = {
                cases[index]["t"]
                for index in indices
            }

            if len(t_values) == 1:
                # Every case in this group is now uniquely
                # identified with respect to t.
                for index in indices:
                    result[index] = k
            else:
                for index in indices:
                    next_unresolved.add(index)

        unresolved = next_unresolved

        if not unresolved:
            break

    return result


def generate_dataset(bits_p, bits_q, count, rng):
    """
    Generate a dataset of distinct semiprimes.
    """
    cases = []
    seen_n = set()

    while len(cases) < count:
        p, q = generate_semiprime(bits_p, bits_q, rng)
        n = p * q

        if n in seen_n:
            continue

        seen_n.add(n)
        cases.append(
            construct_case(p, q)
        )

    return cases


def statistics(values):
    """
    Return min, max, mean and median.
    """
    values = sorted(values)

    if not values:
        return None

    minimum = values[0]
    maximum = values[-1]
    mean = sum(values) / len(values)

    middle = len(values) // 2

    if len(values) % 2 == 0:
        median = (
            values[middle - 1]
            + values[middle]
        ) / 2
    else:
        median = values[middle]

    return minimum, maximum, mean, median


def run_scale(bits_p, bits_q, count, rng):
    """
    Run one scale experiment.
    """
    cases = generate_dataset(
        bits_p,
        bits_q,
        count,
        rng,
    )

    max_bits = max(
        case["bits_N"]
        for case in cases
    )

    k_values = minimal_k_per_case(
        cases,
        max_bits,
    )

    required = [
        k_values[index]
        for index in range(len(cases))
        if index in k_values
    ]

    minimum, maximum, mean, median = statistics(
        required
    )

    fully_deterministic = 0

    for k in range(1, max_bits + 1):
        ok, _ = fingerprint_is_deterministic(
            cases,
            k,
        )

        if ok:
            fully_deterministic = k
            break

    return {
        "cases": cases,
        "N_min_bits": min(
            case["bits_N"]
            for case in cases
        ),
        "N_max_bits": max(
            case["bits_N"]
            for case in cases
        ),
        "min_k": minimum,
        "max_k": maximum,
        "mean_k": mean,
        "median_k": median,
        "full_k": fully_deterministic,
    }


def run_experiment():
    print(f"START EXPERIMENT {EXPERIMENT_NUMBER}")
    print()

    print("2-ADIC SCALING EXPERIMENT")
    print("-------------------------")
    print()
    print("Fingerprint:")
    print("  F_k(N) = (N mod 2^k, floor(sqrt(N)) mod 2^k)")
    print()
    print("Question:")
    print("  Does the amount of 2-adic information needed")
    print("  to determine t remain bounded, or grow with N?")
    print()
    print("For each scale we measure:")
    print("  min k")
    print("  median k")
    print("  mean k")
    print("  max k")
    print("  first k making the entire dataset deterministic")
    print()

    rng = random.Random(32)

    # We use comparable-sized p and q so that the N scale
    # increases in a controlled way.
    #
    # The final entry stays safely below the 64-bit limit.
    scales = [
        (10, 10, 250),
        (14, 14, 250),
        (18, 18, 250),
        (22, 22, 250),
        (26, 26, 250),
        (30, 30, 250),
    ]

    print(
        "bits(p) bits(q) cases  "
        "N-bits        min-k  median-k  mean-k  max-k  full-dataset-k"
    )
    print(
        "----------------------------------------------------------------"
    )

    scale_results = []

    for bits_p, bits_q, count in scales:
        result = run_scale(
            bits_p,
            bits_q,
            count,
            rng,
        )

        scale_results.append(
            (
                bits_p,
                bits_q,
                result,
            )
        )

        print(
            f"{bits_p:7d} "
            f"{bits_q:7d} "
            f"{count:5d}  "
            f"{result['N_min_bits']:3d}.."
            f"{result['N_max_bits']:<3d}    "
            f"{result['min_k']:5d}  "
            f"{result['median_k']:8.1f}  "
            f"{result['mean_k']:7.2f}  "
            f"{result['max_k']:5d}  "
            f"{result['full_k']:14d}"
        )

    print()

    # ---------------------------------------------------------
    # Compare k with log2(N).
    # ---------------------------------------------------------

    print("NORMALIZED k / log2(N)")
    print("----------------------")

    for bits_p, bits_q, result in scale_results:
        cases = result["cases"]

        ratios = []

        # Reconstruct each case's individual required k.
        max_bits = max(
            case["bits_N"]
            for case in cases
        )

        k_values = minimal_k_per_case(
            cases,
            max_bits,
        )

        for index, case in enumerate(cases):
            if index not in k_values:
                continue

            ratios.append(
                k_values[index] / case["bits_N"]
            )

        if ratios:
            print(
                f"{bits_p:2d}+{bits_q:2d} bits: "
                f"min={min(ratios):.4f} "
                f"mean={sum(ratios) / len(ratios):.4f} "
                f"max={max(ratios):.4f}"
            )

    print()

    # ---------------------------------------------------------
    # Test whether k stays below fixed thresholds.
    # ---------------------------------------------------------

    print("FIXED-K COVERAGE")
    print("----------------")

    fixed_k_values = [
        4,
        6,
        8,
        10,
        12,
        14,
        16,
        20,
        24,
    ]

    for fixed_k in fixed_k_values:
        print(
            f"k={fixed_k:2d}:",
            end=" "
        )

        for bits_p, bits_q, result in scale_results:
            cases = result["cases"]

            if fixed_k > max(
                case["bits_N"]
                for case in cases
            ):
                percentage = 100.0
            else:
                k_values = minimal_k_per_case(
                    cases,
                    fixed_k,
                )

                percentage = (
                    100.0
                    * len(k_values)
                    / len(cases)
                )

            print(
                f"{percentage:6.2f}%",
                end=" "
            )

        print()

    print()

    # ---------------------------------------------------------
    # Look for cases requiring unusually large k.
    # ---------------------------------------------------------

    print("LARGEST REQUIRED k CASES")
    print("------------------------")

    all_hard_cases = []

    for bits_p, bits_q, result in scale_results:
        cases = result["cases"]

        max_bits = max(
            case["bits_N"]
            for case in cases
        )

        k_values = minimal_k_per_case(
            cases,
            max_bits,
        )

        for index, k in k_values.items():
            case = cases[index]

            all_hard_cases.append(
                (
                    k,
                    case,
                )
            )

    all_hard_cases.sort(
        key=lambda item: item[0],
        reverse=True
    )

    for k, case in all_hard_cases[:30]:
        print(
            f"N={case['N']} "
            f"p={case['p']} "
            f"q={case['q']} "
            f"Nbits={case['bits_N']} "
            f"t={case['t']} "
            f"required_k={k}"
        )

    print()

    # ---------------------------------------------------------
    # Check for a particularly interesting outcome:
    # Does max required k stay approximately constant?
    # ---------------------------------------------------------

    print("MAX-k SCALING")
    print("-------------")

    previous_max = None

    for bits_p, bits_q, result in scale_results:
        current_max = result["max_k"]

        if previous_max is None:
            delta = "N/A"
        else:
            delta = str(current_max - previous_max)

        print(
            f"N approximately {bits_p + bits_q:2d} bits: "
            f"max_k={current_max:3d} "
            f"delta={delta}"
        )

        previous_max = current_max

    print()

    print("INTERPRETATION")
    print("--------------")
    print("1. The key quantity is max_k as N grows.")
    print("2. If max_k remains roughly constant, that supports")
    print("   a bounded 2-adic fingerprint hypothesis.")
    print("3. If max_k grows with the bit length of N, the earlier")
    print("   2^10 result was probably a scale-dependent effect.")
    print("4. Also compare mean_k and k/log2(N).")
    print("5. The fixed-k coverage table shows whether a practical")
    print("   constant number of low bits survives scaling.")

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT_NUMBER}")


if __name__ == "__main__":
    run_experiment()
