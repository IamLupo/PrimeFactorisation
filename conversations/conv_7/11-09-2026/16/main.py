import math
import random


EXPERIMENT = 38


def is_prime(n):
    if n < 2:
        return False

    small_primes = [
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
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

        composite = True

        for _ in range(r - 1):
            x = (x * x) % n

            if x == n - 1:
                composite = False
                break

        if composite:
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


def valuation_2(n):
    if n == 0:
        return 10**9

    n = abs(n)

    return (n & -n).bit_length() - 1


def polynomial_constant(n, s, d, x):
    # P(x,y) = A(x) + B(x)*y
    #
    # A(x) = x^2 - x + D - s
    # B(x) = x - s - 1

    return (
        x * x
        - x
        + d
        - s
    )


def polynomial_y_coefficient(s, x):
    return x - s - 1


def has_y_solution_mod_power_of_two(n, s, d, x, k):
    """
    Determine whether there exists some y modulo 2^k
    satisfying P(x,y) == 0 mod 2^k.

    P(x,y) = A + B*y

    A + B*y == 0 mod m

    has a solution iff gcd(B,m) divides A.
    """

    modulus = 1 << k

    a = polynomial_constant(n, s, d, x)
    b = polynomial_y_coefficient(s, x)

    g = math.gcd(abs(b), modulus)

    return a % g == 0


def compatible_x_residues(n, s, d, k):
    """
    Return all x mod 4 for which there exists at least
    one x modulo 2^k and some y modulo 2^k satisfying
    P(x,y) == 0 mod 2^k.
    """

    modulus = 1 << k

    residues = set()

    for x in range(modulus):
        if has_y_solution_mod_power_of_two(
            n,
            s,
            d,
            x,
            k,
        ):
            residues.add(x % 4)

            if len(residues) == 4:
                break

    return tuple(sorted(residues))


def exact_two_root_residues(case):
    return tuple(sorted({
        case["x_p"] % 4,
        case["x_q"] % 4,
    }))


def classify_candidate_set(case, candidates):
    actual = exact_two_root_residues(case)

    if not set(actual).issubset(set(candidates)):
        return "MISSING"

    if set(candidates) == set(actual):
        if len(actual) == 1:
            return "EXACT_SINGLE"
        return "EXACT_TWO"

    return "SPURIOUS"


def find_first_example(cases, k, wanted_class):
    for case in cases:
        candidates = compatible_x_residues(
            case["N"],
            case["s"],
            case["D"],
            k,
        )

        classification = classify_candidate_set(
            case,
            candidates,
        )

        if classification == wanted_class:
            return case, candidates

    return None, None


def print_case(case, candidates):
    print(f"N       = {case['N']}")
    print(f"p       = {case['p']}")
    print(f"q       = {case['q']}")
    print(f"s       = {case['s']}")
    print(f"D       = {case['D']}")
    print(f"x_p     = {case['x_p']}")
    print(f"x_q     = {case['x_q']}")
    print(f"x_p%4   = {case['x_p'] % 4}")
    print(f"x_q%4   = {case['x_q'] % 4}")
    print(f"actual  = {exact_two_root_residues(case)}")
    print(f"found   = {candidates}")


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(38038)

    CASE_COUNT = 6000
    PRIME_BITS = 30

    print()
    print("Generating cases...")
    print(f"cases      = {CASE_COUNT}")
    print(f"prime bits = {PRIME_BITS}")
    print()

    cases = []

    for _ in range(CASE_COUNT):
        cases.append(
            generate_case(PRIME_BITS)
        )

    print("============================================================")
    print("TEST 1: P(x,y) ROOT RESIDUES")
    print("============================================================")

    for k in range(2, 13):
        exact_single = 0
        exact_two = 0
        spurious = 0
        missing = 0

        for case in cases:
            candidates = compatible_x_residues(
                case["N"],
                case["s"],
                case["D"],
                k,
            )

            classification = classify_candidate_set(
                case,
                candidates,
            )

            if classification == "EXACT_SINGLE":
                exact_single += 1
            elif classification == "EXACT_TWO":
                exact_two += 1
            elif classification == "SPURIOUS":
                spurious += 1
            elif classification == "MISSING":
                missing += 1

        print(
            f"k={k:2d}  "
            f"single={exact_single:5d}  "
            f"two={exact_two:5d}  "
            f"spurious={spurious:5d}  "
            f"missing={missing:5d}"
        )

    print()

    print("============================================================")
    print("TEST 2: DOES MOD 4 ALREADY GIVE THE TWO TRUE ROOTS?")
    print("============================================================")

    exact = 0
    wrong = 0

    for case in cases:
        candidates = compatible_x_residues(
            case["N"],
            case["s"],
            case["D"],
            2,
        )

        actual = exact_two_root_residues(case)

        if set(candidates) == set(actual):
            exact += 1
        else:
            wrong += 1

    print(
        f"exact candidate set at k=2: "
        f"{exact}/{CASE_COUNT}"
    )

    print(
        f"wrong/spurious candidate set: "
        f"{wrong}/{CASE_COUNT}"
    )

    print()

    print("============================================================")
    print("TEST 3: ROOT SEPARATION")
    print("============================================================")

    same_mod4 = 0
    different_mod4 = 0

    for case in cases:
        if case["x_p"] % 4 == case["x_q"] % 4:
            same_mod4 += 1
        else:
            different_mod4 += 1

    print(
        f"x_p == x_q mod 4: "
        f"{same_mod4}"
    )

    print(
        f"x_p != x_q mod 4: "
        f"{different_mod4}"
    )

    print()

    print("============================================================")
    print("TEST 4: EXAMPLES")
    print("============================================================")

    for wanted in [
        "EXACT_SINGLE",
        "EXACT_TWO",
        "SPURIOUS",
    ]:
        case, candidates = find_first_example(
            cases,
            8,
            wanted,
        )

        print()
        print(f"Example class: {wanted}")

        if case is None:
            print("None found.")
        else:
            print_case(case, candidates)

    print()

    print("============================================================")
    print("TEST 5: DOES MORE 2-ADIC PRECISION REMOVE THE AMBIGUITY?")
    print("============================================================")

    for k in [2, 4, 6, 8, 10, 12]:
        ambiguous = 0
        single_correct = 0

        for case in cases:
            candidates = compatible_x_residues(
                case["N"],
                case["s"],
                case["D"],
                k,
            )

            if len(candidates) > 1:
                ambiguous += 1
            elif (
                len(candidates) == 1
                and candidates[0] == case["x_p"] % 4
            ):
                single_correct += 1

        print(
            f"k={k:2d}  "
            f"ambiguous={ambiguous:5d}  "
            f"unique-correct={single_correct:5d}"
        )

    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()
