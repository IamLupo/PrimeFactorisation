import math

from sympy import isprime


EXPERIMENT = 90


def alternating_transform(x, y, r):
    """
    Exact binomial transform of

        y, x, y, x, ...

    """

    total = 0

    for j in range(r):
        coefficient = (
            math.comb(r - 1, j)
            * (2 ** (r - 1 - j))
        )

        if (r + j) % 2 == 1:
            value = y
        else:
            value = x

        total += coefficient * value

    return total


def alternating_closed_form(x, y, r):
    return (
        (
            (x + y)
            * (3 ** (r - 1))
        )
        +
        (
            (y - x)
            * ((-1) ** (r + 1))
        )
    ) // 2


def verify_closed_form():
    failures = 0

    for x in range(-10, 11):
        for y in range(-10, 11):
            for r in range(1, 15):
                direct = alternating_transform(
                    x,
                    y,
                    r
                )

                closed = alternating_closed_form(
                    x,
                    y,
                    r
                )

                if direct != closed:
                    failures += 1

    return failures


def family_values(n):
    return {
        "n-1,1": (
            n - 1,
            1
        ),

        "n-2,2": (
            n - 2,
            2
        ),

        "n-3,3": (
            n - 3,
            3
        ),

        "n-4,4": (
            n - 4,
            4
        ),

        "n-5,5": (
            n - 5,
            5
        ),

        "n+1,-1": (
            n + 1,
            -1
        ),

        "n+2,-2": (
            n + 2,
            -2
        ),

        "2n-1,1": (
            2 * n - 1,
            1
        ),

        "n,0": (
            n,
            0
        ),

        "1,n-1": (
            1,
            n - 1
        ),
    }


def analyze_family(
    name,
    x,
    y,
    n,
    max_r
):
    print(
        f"{name}: x={x}, y={y}"
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

        gcd_value = math.gcd(
            b,
            n
        )

        gcd_plus = math.gcd(
            b + 1,
            n
        )

        gcd_minus = math.gcd(
            b - 1,
            n
        )

        print(
            f"  r={r:2d} "
            f"b={b} "
            f"b mod n={b % n} "
            f"gcd(b,n)={gcd_value} "
            f"gcd(b-1,n)={gcd_minus} "
            f"gcd(b+1,n)={gcd_plus}"
        )

    print()


def test_special_index(
    x,
    y,
    n
):
    b_n = alternating_closed_form(
        x,
        y,
        n
    )

    return {
        "b_n": b_n,
        "mod_n": b_n % n,
        "gcd": math.gcd(
            b_n,
            n
        ),
        "gcd_minus": math.gcd(
            b_n - 1,
            n
        ),
        "gcd_plus": math.gcd(
            b_n + 1,
            n
        ),
        "prime_b_n": isprime(b_n)
    }


def scan_prime_composite(
    max_n=200
):
    families = [
        (
            "n-1,1",
            lambda n: (
                n - 1,
                1
            )
        ),
        (
            "n-2,2",
            lambda n: (
                n - 2,
                2
            )
        ),
        (
            "n-3,3",
            lambda n: (
                n - 3,
                3
            )
        ),
        (
            "n-4,4",
            lambda n: (
                n - 4,
                4
            )
        ),
        (
            "n+1,-1",
            lambda n: (
                n + 1,
                -1
            )
        ),
        (
            "2n-1,1",
            lambda n: (
                2 * n - 1,
                1
            )
        ),
    ]

    for name, function in families:
        prime_hits = 0
        prime_misses = 0
        composite_factor_hits = 0

        for n in range(
            2,
            max_n + 1
        ):
            x, y = function(n)

            result = test_special_index(
                x,
                y,
                n
            )

            if isprime(n):
                if result["gcd"] == 1:
                    prime_hits += 1
                else:
                    prime_misses += 1
            else:
                if (
                    result["gcd"] != 1
                    and result["gcd"] != n
                ):
                    composite_factor_hits += 1

        print(name)
        print(
            f"  PRIME n, gcd(b_n,n)=1: "
            f"{prime_hits}"
        )
        print(
            f"  PRIME n, gcd(b_n,n)!=1: "
            f"{prime_misses}"
        )
        print(
            f"  COMPOSITE n, nontrivial gcd: "
            f"{composite_factor_hits}"
        )
        print()


def search_linear_parameters(
    max_abs=5,
    max_n=100
):
    """
    Search

        x = a*n + b
        y = c*n + d

    for small coefficients.

    We specifically look for families where

        gcd(b_n, n)

    behaves differently for primes and composites.
    """

    candidates = []

    for a in range(
        -max_abs,
        max_abs + 1
    ):
        for b in range(
            -max_abs,
            max_abs + 1
        ):
            for c in range(
                -max_abs,
                max_abs + 1
            ):
                for d in range(
                    -max_abs,
                    max_abs + 1
                ):

                    prime_good = 0
                    composite_factor = 0

                    for n in range(
                        3,
                        max_n + 1
                    ):
                        x = a * n + b
                        y = c * n + d

                        value = alternating_closed_form(
                            x,
                            y,
                            n
                        )

                        g = math.gcd(
                            value,
                            n
                        )

                        if isprime(n):
                            if g == 1:
                                prime_good += 1

                        else:
                            if (
                                g != 1
                                and g != n
                            ):
                                composite_factor += 1

                    if prime_good >= 20:
                        candidates.append(
                            (
                                composite_factor,
                                prime_good,
                                a,
                                b,
                                c,
                                d
                            )
                        )

    candidates.sort(
        reverse=True
    )

    return candidates[:20]


def run_experiment():
    print("=" * 60)
    print(
        f"START EXPERIMENT {EXPERIMENT}"
    )
    print("=" * 60)
    print()

    print(
        "1. CLOSED FORM VERIFICATION"
    )
    print("-" * 60)

    failures = verify_closed_form()

    print(
        f"FAILURES: {failures}"
    )
    print()

    print(
        "2. EXPLICIT ALTERNATING FAMILIES"
    )
    print("-" * 60)

    n = 77

    families = family_values(n)

    for name, (
        x,
        y
    ) in families.items():
        analyze_family(
            name,
            x,
            y,
            n,
            max_r=8
        )

    print(
        "3. PRIME/COMPOSITE SCAN"
    )
    print("-" * 60)

    scan_prime_composite(
        max_n=200
    )

    print(
        "4. LINEAR FAMILY SEARCH"
    )
    print("-" * 60)

    candidates = search_linear_parameters(
        max_abs=3,
        max_n=100
    )

    if not candidates:
        print(
            "NO CANDIDATES FOUND"
        )
    else:
        print(
            "TOP CANDIDATES:"
        )

        for candidate in candidates:
            (
                composite_factor,
                prime_good,
                a,
                b,
                c,
                d
            ) = candidate

            print(
                "  "
                f"x={a}n+{b}, "
                f"y={c}n+{d}, "
                f"prime_good={prime_good}, "
                f"composite_factor={composite_factor}"
            )

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
