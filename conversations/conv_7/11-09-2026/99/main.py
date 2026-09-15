import math
import random
import time

from sympy import isprime


EXPERIMENT = 122

CASES = 50

P_MIN = 100_000
P_MAX = 1_000_000

Q_MIN = 1_000_000
Q_MAX = 10_000_000

MAX_PROBES = 16

# Compare exact math.comb() against the modular evaluator
# on only this many cases. After that, only the fast
# evaluator is used.
VALIDATION_CASES = 10


def gcd_state(g, n, p, q):
    if g == 1:
        return "1"

    if g == n:
        return "n"

    if g == p:
        return "p"

    if g == q:
        return "q"

    return "other"


def bit_reverse(value, bits):
    result = 0

    for _ in range(bits):
        result <<= 1
        result |= value & 1
        value >>= 1

    return result


def van_der_corput_sequence(max_z, count):
    count = min(count, max_z)

    if count <= 0:
        return []

    bits = max(
        1,
        (count - 1).bit_length(),
    )

    result = []
    used = set()

    index = 0

    while len(result) < count:
        reversed_value = bit_reverse(
            index,
            bits,
        )

        z = (
            reversed_value * max_z
            // (1 << bits)
        ) + 1

        z = max(
            1,
            min(z, max_z),
        )

        if z not in used:
            used.add(z)
            result.append(z)

        index += 1

    return result


def exact_binomial_mod(n, k, z):
    """
    Reference implementation.

    Constructs the full binomial coefficient.
    """

    d = 2 * z - 1

    if d >= k - 1:
        return n, None, 0.0

    start = time.perf_counter()

    value = math.comb(
        k - 2,
        d,
    )

    value %= n
    value = (-value) % n

    elapsed = (
        time.perf_counter()
        - start
    )

    return value, None, elapsed


def modular_binomial_mod(n, k, z):
    """
    Compute

        C(k-2, 2z-1) mod n

    using modular inverses.

    If some denominator j is not invertible modulo n,
    then gcd(j,n) gives a nontrivial factor immediately.

    Returns:

        value_mod_n
        factor_or_none
        elapsed
        iterations
    """

    d = 2 * z - 1

    if d >= k - 1:
        return n, None, 0.0, 0

    m = k - 2

    start = time.perf_counter()

    result = 1

    for j in range(1, d + 1):
        g = math.gcd(
            j,
            n,
        )

        if 1 < g < n:
            elapsed = (
                time.perf_counter()
                - start
            )

            return (
                None,
                g,
                elapsed,
                j,
            )

        if g != 1:
            # gcd(j,n) == n cannot happen here because
            # j < n in all intended test cases.
            return (
                None,
                g,
                time.perf_counter() - start,
                j,
            )

        numerator = m - d + j

        result *= numerator % n
        result %= n

        inverse = pow(
            j,
            -1,
            n,
        )

        result *= inverse
        result %= n

    elapsed = (
        time.perf_counter()
        - start
    )

    return (
        result,
        None,
        elapsed,
        d,
    )


def value_gcd_exact(n, k, z):
    value, factor, elapsed = (
        exact_binomial_mod(
            n,
            k,
            z,
        )
    )

    if factor is not None:
        return (
            factor,
            elapsed,
        )

    value = (-value) % n

    return (
        math.gcd(
            value,
            n,
        ),
        elapsed,
    )


def value_gcd_modular(n, k, z):
    (
        value,
        factor,
        elapsed,
        iterations,
    ) = modular_binomial_mod(
        n,
        k,
        z,
    )

    if factor is not None:
        return (
            factor,
            elapsed,
            iterations,
        )

    value = (-value) % n

    return (
        math.gcd(
            value,
            n,
        ),
        elapsed,
        iterations,
    )


def generate_semiprime():
    while True:
        p = random.randrange(
            P_MIN,
            P_MAX + 1,
        )

        if not isprime(p):
            continue

        q = random.randrange(
            max(p + 1, Q_MIN),
            Q_MAX + 1,
        )

        if not isprime(q):
            continue

        return p, q


def powers_between(p, q):
    result = []

    k = 1

    while k <= p:
        k *= 2

    while k < q:
        result.append(k)
        k *= 2

    return result


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {CASES}")
    print(f"P RANGE: {P_MIN}..{P_MAX}")
    print(f"Q RANGE: {Q_MIN}..{Q_MAX}")
    print(f"MAX PROBES: {MAX_PROBES}")
    print(f"VALIDATION CASES: {VALIDATION_CASES}")
    print()

    random.seed(EXPERIMENT)

    exact_total_time = 0.0
    modular_total_time = 0.0

    exact_evaluations = 0
    modular_evaluations = 0

    modular_factors_found = 0

    validation_matches = 0
    validation_mismatches = 0

    total_cases = 0

    exact_hits = 0
    modular_hits = 0

    modular_iterations = 0

    max_exact_time = 0.0
    max_modular_time = 0.0

    max_exact_example = None
    max_modular_example = None

    mismatch_examples = []

    for case_index in range(
        1,
        CASES + 1,
    ):
        p, q = generate_semiprime()
        n = p * q

        ks = powers_between(
            p,
            q,
        )

        if not ks:
            continue

        # Keep the same first useful k as Experiment 121.
        k = ks[0]

        max_z = k // 2

        probes = van_der_corput_sequence(
            max_z,
            MAX_PROBES,
        )

        total_cases += 1

        exact_found = False
        modular_found = False

        for probe_index, z in enumerate(
            probes,
            start=1,
        ):
            # -------------------------------------------------
            # MODULAR
            # -------------------------------------------------

            (
                modular_g,
                modular_time,
                iterations,
            ) = value_gcd_modular(
                n,
                k,
                z,
            )

            modular_total_time += (
                modular_time
            )

            modular_evaluations += 1
            modular_iterations += (
                iterations
            )

            if modular_time > max_modular_time:
                max_modular_time = (
                    modular_time
                )

                max_modular_example = (
                    case_index,
                    n,
                    p,
                    q,
                    k,
                    z,
                    iterations,
                )

            if (
                modular_g is not None
                and 1 < modular_g < n
            ):
                modular_hits += 1
                modular_found = True
                modular_factors_found += 1

            # -------------------------------------------------
            # EXACT REFERENCE
            # -------------------------------------------------

            if (
                case_index
                <= VALIDATION_CASES
            ):
                (
                    exact_g,
                    exact_time,
                ) = value_gcd_exact(
                    n,
                    k,
                    z,
                )

                exact_total_time += (
                    exact_time
                )

                exact_evaluations += 1

                if exact_time > max_exact_time:
                    max_exact_time = (
                        exact_time
                    )

                    max_exact_example = (
                        case_index,
                        n,
                        p,
                        q,
                        k,
                        z,
                    )

                if modular_g != exact_g:
                    validation_mismatches += 1

                    if len(
                        mismatch_examples
                    ) < 10:
                        mismatch_examples.append(
                            (
                                case_index,
                                n,
                                p,
                                q,
                                k,
                                z,
                                exact_g,
                                modular_g,
                            )
                        )
                else:
                    validation_matches += 1

            # Stop once this strategy found a factor.
            if modular_found:
                break

        if exact_found:
            exact_hits += 1

    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(
        f"cases                       = "
        f"{total_cases}"
    )

    print(
        f"modular factor hits         = "
        f"{modular_factors_found}"
    )

    print(
        f"validation matches          = "
        f"{validation_matches}"
    )

    print(
        f"validation mismatches       = "
        f"{validation_mismatches}"
    )

    print()

    print(
        f"modular total time          = "
        f"{modular_total_time:.6f}s"
    )

    print(
        f"modular average evaluation = "
        f"{modular_total_time / max(1, modular_evaluations):.6f}s"
    )

    if exact_evaluations:
        print(
            f"exact total time            = "
            f"{exact_total_time:.6f}s"
        )

        print(
            f"exact average evaluation   = "
            f"{exact_total_time / exact_evaluations:.6f}s"
        )

        print(
            f"speedup                     = "
            f"{exact_total_time / max(modular_total_time, 1e-12):.2f}x"
        )

    print()

    print(
        f"modular iterations          = "
        f"{modular_iterations}"
    )

    if modular_evaluations:
        print(
            f"average modular iterations = "
            f"{modular_iterations / modular_evaluations:.2f}"
        )

    print()

    print(
        f"maximum exact evaluation    = "
        f"{max_exact_time:.6f}s"
    )

    print(
        f"maximum modular evaluation  = "
        f"{max_modular_time:.6f}s"
    )

    if max_exact_example is not None:
        (
            case_index,
            n,
            p,
            q,
            k,
            z,
        ) = max_exact_example

        print()
        print("SLOWEST EXACT EVALUATION")

        print(
            f"  CASE={case_index}"
        )

        print(
            f"  N={n}"
        )

        print(
            f"  p={p}"
        )

        print(
            f"  q={q}"
        )

        print(
            f"  k={k}"
        )

        print(
            f"  z={z}"
        )

    if max_modular_example is not None:
        (
            case_index,
            n,
            p,
            q,
            k,
            z,
            iterations,
        ) = max_modular_example

        print()
        print("SLOWEST MODULAR EVALUATION")

        print(
            f"  CASE={case_index}"
        )

        print(
            f"  N={n}"
        )

        print(
            f"  p={p}"
        )

        print(
            f"  q={q}"
        )

        print(
            f"  k={k}"
        )

        print(
            f"  z={z}"
        )

        print(
            f"  iterations={iterations}"
        )

    print()
    print("-" * 60)
    print("VALIDATION")
    print("-" * 60)

    if validation_mismatches == 0:
        print("ALL MODULAR RESULTS MATCH EXACT RESULTS.")
    else:
        print(
            f"MISMATCHES={validation_mismatches}"
        )

        for example in mismatch_examples:
            print(example)

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
