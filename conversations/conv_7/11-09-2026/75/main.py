import math
import random


EXPERIMENT = 94


def is_prime(n):
    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    limit = math.isqrt(n)

    d = 3

    while d <= limit:
        if n % d == 0:
            return False

        d += 2

    return True


def generate_semiprime(
    min_prime=11,
    max_prime=97
):
    while True:
        p = random.randrange(
            min_prime,
            max_prime + 1
        )

        q = random.randrange(
            min_prime,
            max_prime + 1
        )

        if p == q:
            continue

        if not is_prime(p):
            continue

        if not is_prime(q):
            continue

        if p > q:
            p, q = q, p

        return p * q, p, q


def arthseq_value(
    n,
    repetition,
    k
):
    degree = 2 * repetition - 1

    total = 0

    for j in range(
        degree + 1
    ):
        if j % 2 == 0:
            value = 1
        else:
            value = n - 1

        total += (
            value
            * math.comb(
                k - 1,
                j
            )
        )

    return total


def predicted_prime_hits(
    prime,
    degree,
    count
):
    hits = set()

    block_width = degree

    multiple = 1

    while True:
        start = (
            multiple * prime
            + 2
        )

        if start > count:
            break

        for offset in range(
            block_width
        ):
            position = (
                start
                + offset
            )

            if position > count:
                break

            hits.add(position)

        multiple += 1

    return hits


def actual_prime_hits(
    n,
    prime,
    repetition,
    count
):
    hits = set()

    for k in range(
        1,
        count + 1
    ):
        value = arthseq_value(
            n,
            repetition,
            k
        )

        if value % prime == 0:
            hits.add(k)

    return hits


def verify_case(
    n,
    p,
    q,
    repetition,
    count
):
    degree = (
        2 * repetition - 1
    )

    predicted_p = predicted_prime_hits(
        p,
        degree,
        count
    )

    predicted_q = predicted_prime_hits(
        q,
        degree,
        count
    )

    actual_p = actual_prime_hits(
        n,
        p,
        repetition,
        count
    )

    actual_q = actual_prime_hits(
        n,
        q,
        repetition,
        count
    )

    return (
        predicted_p,
        actual_p,
        predicted_q,
        actual_q
    )


def run_experiment(
    num_cases=500,
    max_repetition=10,
    count=100
):
    print("=" * 60)
    print(
        f"START EXPERIMENT {EXPERIMENT}"
    )
    print("=" * 60)
    print()

    print(
        f"CASES: {num_cases}"
    )

    print(
        f"FIRST VALUES: {count}"
    )

    print()

    total_p_failures = 0
    total_q_failures = 0

    total_valid_tests = 0

    block_width_failures = 0

    printed_cases = 0

    for case in range(
        1,
        num_cases + 1
    ):
        n, p, q = generate_semiprime()

        for repetition in range(
            1,
            max_repetition + 1
        ):
            degree = (
                2 * repetition - 1
            )

            # Stay strictly inside the regime
            # degree < p and degree < q.
            if degree >= p:
                continue

            if degree >= q:
                continue

            total_valid_tests += 1

            (
                predicted_p,
                actual_p,
                predicted_q,
                actual_q
            ) = verify_case(
                n,
                p,
                q,
                repetition,
                count
            )

            if predicted_p != actual_p:
                total_p_failures += 1

            if predicted_q != actual_q:
                total_q_failures += 1

            # Verify that each contiguous block
            # has exactly degree positions unless
            # it is cut off by the 100-value window.
            for prime in (
                p,
                q
            ):
                predicted = predicted_prime_hits(
                    prime,
                    degree,
                    count
                )

                actual = actual_prime_hits(
                    n,
                    prime,
                    repetition,
                    count
                )

                if predicted != actual:
                    block_width_failures += 1

        if printed_cases < 10:
            printed_cases += 1

            print(
                f"CASE {printed_cases}"
            )

            print(
                f"  N = {n}"
            )

            print(
                f"  p = {p}"
            )

            print(
                f"  q = {q}"
            )

            for repetition in range(
                1,
                max_repetition + 1
            ):
                degree = (
                    2 * repetition - 1
                )

                if degree >= p:
                    break

                (
                    predicted_p,
                    actual_p,
                    predicted_q,
                    actual_q
                ) = verify_case(
                    n,
                    p,
                    q,
                    repetition,
                    count
                )

                print(
                    f"  L={repetition:2d} "
                    f"d={degree:2d} "
                    f"p_match="
                    f"{predicted_p == actual_p} "
                    f"q_match="
                    f"{predicted_q == actual_q}"
                )

                print(
                    f"      p="
                    f"{sorted(actual_p)}"
                )

                print(
                    f"      q="
                    f"{sorted(actual_q)}"
                )

            print()

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)

    print(
        f"VALID TESTS: "
        f"{total_valid_tests}"
    )

    print(
        f"p PREDICTION FAILURES: "
        f"{total_p_failures}"
    )

    print(
        f"q PREDICTION FAILURES: "
        f"{total_q_failures}"
    )

    print(
        f"BLOCK WIDTH FAILURES: "
        f"{block_width_failures}"
    )

    print()

    if (
        total_p_failures == 0
        and total_q_failures == 0
    ):
        print(
            "RESULT: EXACT BLOCK LAW VERIFIED"
        )

        print(
            "For degree d=2L-1 < p,q:"
        )

        print(
            "  v_k == 0 (mod p)"
            " exactly at predicted blocks."
        )

        print(
            "  v_k == 0 (mod q)"
            " exactly at predicted blocks."
        )

    print()
    print("=" * 60)
    print(
        f"FINISHED EXPERIMENT {EXPERIMENT}"
    )
    print("=" * 60)


def main():
    random.seed(940011)

    run_experiment(
        num_cases=500,
        max_repetition=10,
        count=100
    )


if __name__ == "__main__":
    main()
