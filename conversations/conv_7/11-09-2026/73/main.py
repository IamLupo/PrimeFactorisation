import math

from sympy import isprime


EXPERIMENT = 91


def alternating_closed_form(x, y, r):
    numerator = (
        (x + y) * (3 ** (r - 1))
        +
        (y - x) * ((-1) ** (r + 1))
    )

    return numerator // 2


def family_value(n, c):
    x = n - c
    y = c

    return x, y


def verify_family_congruence(
    max_n=200,
    max_c=20
):
    failures = 0

    for n in range(
        3,
        max_n + 1
    ):
        if n % 2 == 0:
            continue

        for c in range(
            -max_c,
            max_c + 1
        ):
            x, y = family_value(
                n,
                c
            )

            b = alternating_closed_form(
                x,
                y,
                n
            )

            expected = (
                c
                * ((-1) ** (n + 1))
            ) % n

            if b % n != expected:
                failures += 1

    return failures


def scan_indices(
    n,
    c,
    max_r
):
    x, y = family_value(
        n,
        c
    )

    result = []

    for r in range(
        1,
        max_r + 1
    ):
        b = alternating_closed_form(
            x,
            y,
            r
        )

        g = math.gcd(
            b,
            n
        )

        if (
            g != 1
            and g != n
        ):
            result.append(
                (
                    r,
                    b,
                    g
                )
            )

    return result


def semiprime(p, q):
    return p * q


def generate_semiprimes(
    min_prime=100,
    max_prime=5000,
    count=200
):
    values = []

    p = min_prime

    while p <= max_prime:
        if isprime(p):
            q = p + 1

            while q <= max_prime:
                if (
                    isprime(q)
                    and p != q
                ):
                    values.append(
                        (
                            p * q,
                            p,
                            q
                        )
                    )

                    if len(values) >= count:
                        return values

                q += 1

        p += 1

    return values


def measure_family(
    c,
    semiprimes,
    max_r
):
    cases_with_factor = 0
    total_hits = 0

    for n, p, q in semiprimes:
        hits = scan_indices(
            n,
            c,
            max_r
        )

        if hits:
            cases_with_factor += 1
            total_hits += len(hits)

    return (
        cases_with_factor,
        total_hits
    )


def compare_c_values():
    semiprimes = generate_semiprimes(
        min_prime=100,
        max_prime=1000,
        count=200
    )

    print(
        f"SEMIPRIME CASES: "
        f"{len(semiprimes)}"
    )

    for c in range(
        0,
        11
    ):
        (
            cases_with_factor,
            total_hits
        ) = measure_family(
            c,
            semiprimes,
            max_r=64
        )

        print(
            f"c={c:2d} "
            f"factor cases="
            f"{cases_with_factor}/"
            f"{len(semiprimes)} "
            f"total hits="
            f"{total_hits}"
        )


def show_case(
    n,
    c,
    max_r=32
):
    x, y = family_value(
        n,
        c
    )

    print(
        f"N={n} "
        f"c={c} "
        f"x={x} "
        f"y={y}"
    )

    for r in range(
        1,
        max_r + 1
    ):
        b = alternating_closed_form(
            x,
            y,
            r
        )

        g = math.gcd(
            b,
            n
        )

        if (
            g != 1
            and g != n
        ):
            print(
                f"  r={r} "
                f"gcd={g} "
                f"b={b}"
            )


def run_experiment():
    print("=" * 60)
    print(
        f"START EXPERIMENT {EXPERIMENT}"
    )
    print("=" * 60)
    print()

    print(
        "1. VERIFY b_n CONGRUENCE"
    )
    print("-" * 60)

    failures = verify_family_congruence()

    print(
        f"FAILURES: {failures}"
    )
    print()

    print(
        "2. EXPLICIT TEST"
    )
    print("-" * 60)

    for c in range(
        0,
        6
    ):
        n = 77

        x, y = family_value(
            n,
            c
        )

        b = alternating_closed_form(
            x,
            y,
            n
        )

        print(
            f"c={c} "
            f"x={x} "
            f"y={y} "
            f"b_n mod n={b % n} "
            f"gcd={math.gcd(b,n)}"
        )

    print()

    print(
        "3. SEARCH OTHER r"
    )
    print("-" * 60)

    for c in range(
        0,
        6
    ):
        print(
            f"c={c}"
        )

        show_case(
            77,
            c,
            max_r=32
        )

        print()

    print(
        "4. SEMIPRIME SCAN"
    )
    print("-" * 60)

    compare_c_values()

    print()

    print("=" * 60)
    print(
        f"FINISHED EXPERIMENT {EXPERIMENT}"
    )
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
