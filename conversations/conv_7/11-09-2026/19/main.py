import math
import random


EXPERIMENT = 40


def is_prime(n):
    if n < 2:
        return False

    small_primes = [
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37
    ]

    for p in small_primes:
        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    r = 0

    while d % 2 == 0:
        d //= 2
        r += 1

    witnesses = [2, 3, 5, 7, 11, 13, 17]

    for a in witnesses:
        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(r - 1):
            x = (x * x) % n

            if x == n - 1:
                break
        else:
            return False

    return True


def random_prime(bits):
    while True:
        n = random.getrandbits(bits)

        n |= 1
        n |= 1 << (bits - 1)

        if is_prime(n):
            return n


def generate_case(bits):
    while True:
        p = random_prime(bits)
        q = random_prime(bits)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q
        s = math.isqrt(n)

        if s <= 4096:
            continue

        d = n - s * s

        x_p = s + 1 - p
        x_q = s + 1 - q

        y = p + q - 2 * s - 1

        return {
            "p": p,
            "q": q,
            "N": n,
            "s": s,
            "D": d,
            "x_p": x_p,
            "x_q": x_q,
            "y": y,
        }


def gcd_condition(n, s, d, x, k):
    modulus = 1 << k

    a = (
        x * x
        - x
        + d
        - s
    )

    b = x - s - 1

    g = math.gcd(abs(b), modulus)

    return a % g == 0


def theoretical_condition(n, s, x):
    """
    Since N is odd:

        gcd(x-s-1, 2^k) | P_constant

    iff x-s-1 is odd.

    Therefore x == s (mod 2).
    """

    return (x - s - 1) % 2 == 1


def actual_x_candidate_count(s, k):
    """
    Exactly half of all residues modulo 2^k satisfy
        x == s (mod 2).
    """

    return 1 << (k - 1)


def positive_candidate_count(s, k):
    """
    Number of positive integers x in [1,s] with
        x == s (mod 2).
    """

    if s <= 0:
        return 0

    if s % 2 == 0:
        return s // 2
    else:
        return (s + 1) // 2


def centered_residue_count(s, k):
    """
    Number of centered representatives satisfying
        x == s (mod 2).
    """

    modulus = 1 << k

    # Exactly half of the residue classes.
    return modulus // 2


def verify_condition(cases):
    print("============================================================")
    print("TEST 1: EXACT EXISTENCE CONDITION")
    print("============================================================")

    for k in range(2, 17):
        failures = 0
        tested = 0

        modulus = 1 << k

        # We only test a small number of residues per case.
        # The algebra says the result must depend only on parity.
        sample_x = []

        for x in range(min(modulus, 32)):
            sample_x.append(x)

        for case in cases[:100]:
            n = case["N"]
            s = case["s"]
            d = case["D"]

            for x in sample_x:
                actual = gcd_condition(
                    n,
                    s,
                    d,
                    x,
                    k,
                )

                predicted = theoretical_condition(
                    n,
                    s,
                    x,
                )

                if actual != predicted:
                    failures += 1

                tested += 1

        print(
            f"k={k:2d}  "
            f"tested={tested:6d}  "
            f"failures={failures:4d}"
        )

    print()


def test_true_roots(cases):
    print("============================================================")
    print("TEST 2: TRUE FACTOR ROOTS")
    print("============================================================")

    for k in range(2, 17):
        xp_ok = 0
        xq_ok = 0

        for case in cases:
            n = case["N"]
            s = case["s"]
            d = case["D"]

            xp = case["x_p"]
            xq = case["x_q"]

            if gcd_condition(n, s, d, xp, k):
                xp_ok += 1

            if gcd_condition(n, s, d, xq, k):
                xq_ok += 1

        print(
            f"k={k:2d}  "
            f"x_p valid={xp_ok:5d}/{len(cases)}  "
            f"x_q valid={xq_ok:5d}/{len(cases)}"
        )

    print()


def test_candidate_counts(cases):
    print("============================================================")
    print("TEST 3: NUMBER OF 2-ADIC X CANDIDATES")
    print("============================================================")

    for k in range(2, 17):
        modulus = 1 << k

        theoretical = modulus // 2

        # All x residues satisfying x == s (mod 2).
        #
        # This is independent of N apart from the parity of s.
        observed = []

        for case in cases[:100]:
            s = case["s"]

            count = sum(
                1
                for x in range(modulus)
                if (x - s) % 2 == 0
            )

            observed.append(count)

        minimum = min(observed)
        maximum = max(observed)

        print(
            f"k={k:2d}  "
            f"theoretical={theoretical:6d}  "
            f"observed={minimum:6d}..{maximum:6d}"
        )

    print()


def test_positive_region(cases):
    print("============================================================")
    print("TEST 4: POSITIVE CANDIDATES")
    print("============================================================")

    for k in [4, 8, 12, 16]:
        counts = []

        for case in cases:
            s = case["s"]

            counts.append(
                positive_candidate_count(
                    s,
                    k,
                )
            )

        print(
            f"k={k:2d}  "
            f"min={min(counts):8d}  "
            f"max={max(counts):8d}  "
            f"mean={sum(counts) / len(counts):.2f}"
        )

    print()


def compare_true_root_parity(cases):
    print("============================================================")
    print("TEST 5: TRUE ROOT PARITY")
    print("============================================================")

    xp_match = 0
    xq_match = 0

    for case in cases:
        s = case["s"]

        xp = case["x_p"]
        xq = case["x_q"]

        if xp % 2 == s % 2:
            xp_match += 1

        if xq % 2 == s % 2:
            xq_match += 1

    print(
        f"x_p == s mod 2: "
        f"{xp_match}/{len(cases)}"
    )

    print(
        f"x_q == s mod 2: "
        f"{xq_match}/{len(cases)}"
    )

    print()


def print_example(case):
    s = case["s"]

    print(f"N     = {case['N']}")
    print(f"p     = {case['p']}")
    print(f"q     = {case['q']}")
    print(f"s     = {s}")
    print(f"x_p   = {case['x_p']}")
    print(f"x_q   = {case['x_q']}")
    print(f"x_p%2 = {case['x_p'] % 2}")
    print(f"x_q%2 = {case['x_q'] % 2}")
    print(f"s%2   = {s % 2}")


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(40040)

    CASE_COUNT = 1500
    PRIME_BITS = 20

    print()
    print("Generating cases...")
    print(f"cases      = {CASE_COUNT}")
    print(f"prime bits = {PRIME_BITS}")
    print()

    cases = [
        generate_case(PRIME_BITS)
        for _ in range(CASE_COUNT)
    ]

    verify_condition(cases)
    test_true_roots(cases)
    test_candidate_counts(cases)
    test_positive_region(cases)
    compare_true_root_parity(cases)

    print("============================================================")
    print("EXAMPLE")
    print("============================================================")

    print_example(cases[0])

    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
