import math
import random


EXPERIMENT = 41


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

    for a in [2, 3, 5, 7, 11, 13, 17]:
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


def polynomial(n, s, d, x, y):
    return (
        x * x
        + (y - 1) * x
        + (d - s)
        - (s + 1) * y
    )


def centered_residue(value, modulus):
    value %= modulus

    half = modulus // 2

    if value >= half:
        value -= modulus

    return value


def solve_y_mod(n, s, d, x, modulus):
    """
    P(x,y) = A + B*y

    A = x^2 - x + D - s
    B = x - s - 1

    For admissible x, B is odd and therefore invertible modulo 2^k.
    """

    a = (
        x * x
        - x
        + d
        - s
    )

    b = x - s - 1

    if (b & 1) == 0:
        raise ValueError(
            "solve_y_mod called with even denominator"
        )

    y = (
        -a * pow(b, -1, modulus)
    ) % modulus

    return centered_residue(y, modulus)


def enumerate_candidates(case, k):
    """
    Enumerate all centered x residues modulo 2^k satisfying
    the necessary parity condition x == s (mod 2).

    For every such x there is a unique y modulo 2^k.
    """

    modulus = 1 << k
    half = modulus // 2

    candidates = []

    for x in range(-half, half):
        if ((x - case["s"]) & 1) != 0:
            continue

        y = solve_y_mod(
            case["N"],
            case["s"],
            case["D"],
            x,
            modulus,
        )

        residue = polynomial(
            case["N"],
            case["s"],
            case["D"],
            x,
            y,
        )

        candidates.append({
            "x": x,
            "y": y,
            "residue": residue,
            "residue_mod": residue % modulus,
        })

    return candidates


def geometric_candidates(candidates, s):
    """
    Candidate geometry based on the transformed factor coordinates:

        p = s + 1 - x
        q = s + x + y

    We require:

        p >= 2
        q > s
        y >= 1
    """

    result = []

    for candidate in candidates:
        x = candidate["x"]
        y = candidate["y"]

        p_candidate = s + 1 - x
        q_candidate = s + x + y

        if (
            p_candidate >= 2
            and q_candidate > s
            and y >= 1
        ):
            result.append({
                **candidate,
                "p_candidate": p_candidate,
                "q_candidate": q_candidate,
            })

    return result


def true_residue_pair(case, k):
    modulus = 1 << k

    return (
        case["x_p"] % modulus,
        case["y"] % modulus,
    )


def centered_true_residue_pair(case, k):
    modulus = 1 << k

    return (
        centered_residue(
            case["x_p"],
            modulus,
        ),
        centered_residue(
            case["y"],
            modulus,
        ),
    )


def contains_true_residue_pair(candidates, case, k):
    true_x, true_y = centered_true_residue_pair(
        case,
        k,
    )

    for candidate in candidates:
        if (
            candidate["x"] == true_x
            and candidate["y"] == true_y
        ):
            return True

    return False


def contains_true_x_residue(candidates, case, k):
    true_x = centered_residue(
        case["x_p"],
        1 << k,
    )

    return any(
        candidate["x"] == true_x
        for candidate in candidates
    )


def geometric_contains_true_residue(candidates, case, k):
    true_x, true_y = centered_true_residue_pair(
        case,
        k,
    )

    return any(
        candidate["x"] == true_x
        and candidate["y"] == true_y
        for candidate in candidates
    )


def exact_zero(candidates):
    return [
        candidate
        for candidate in candidates
        if candidate["residue"] == 0
    ]


def print_case(case, k):
    modulus = 1 << k

    true_x = centered_residue(
        case["x_p"],
        modulus,
    )

    true_y = centered_residue(
        case["y"],
        modulus,
    )

    print(f"N             = {case['N']}")
    print(f"p             = {case['p']}")
    print(f"q             = {case['q']}")
    print(f"s             = {case['s']}")
    print(f"D             = {case['D']}")
    print(f"x_p           = {case['x_p']}")
    print(f"y_true        = {case['y']}")
    print(f"modulus       = {modulus}")
    print(f"x_p mod M     = {case['x_p'] % modulus}")
    print(f"y_true mod M   = {case['y'] % modulus}")
    print(f"centered x    = {true_x}")
    print(f"centered y    = {true_y}")


def main():
    print(f"START EXPERIMENT {EXPERIMENT}")

    random.seed(414141)

    CASE_COUNT = 100
    PRIME_BITS = 14

    print()
    print("Generating cases...")
    print(f"cases      = {CASE_COUNT}")
    print(f"prime bits = {PRIME_BITS}")
    print()

    cases = [
        generate_case(PRIME_BITS)
        for _ in range(CASE_COUNT)
    ]

    for k in [8, 10, 12, 14, 16]:
        print("============================================================")
        print(f"TEST k={k}")
        print("============================================================")

        total_candidates = 0
        total_geometric = 0
        total_zero = 0

        true_pair_present = 0
        true_x_present = 0
        true_pair_geometric = 0

        unique_geometric = 0
        unique_zero = 0

        for case in cases:
            candidates = enumerate_candidates(
                case,
                k,
            )

            geometric = geometric_candidates(
                candidates,
                case["s"],
            )

            zeros = exact_zero(candidates)

            total_candidates += len(candidates)
            total_geometric += len(geometric)
            total_zero += len(zeros)

            if contains_true_residue_pair(
                candidates,
                case,
                k,
            ):
                true_pair_present += 1

            if contains_true_x_residue(
                candidates,
                case,
                k,
            ):
                true_x_present += 1

            if geometric_contains_true_residue(
                geometric,
                case,
                k,
            ):
                true_pair_geometric += 1

            if len(geometric) == 1:
                unique_geometric += 1

            if len(zeros) == 1:
                unique_zero += 1

        print(
            f"total modular candidates   = "
            f"{total_candidates}"
        )

        print(
            f"total geometric candidates = "
            f"{total_geometric}"
        )

        print(
            f"total exact integer zeros   = "
            f"{total_zero}"
        )

        print(
            f"true (x,y) residue present  = "
            f"{true_pair_present}/{CASE_COUNT}"
        )

        print(
            f"true x residue present      = "
            f"{true_x_present}/{CASE_COUNT}"
        )

        print(
            f"true pair survives geometry = "
            f"{true_pair_geometric}/{CASE_COUNT}"
        )

        print(
            f"unique geometric            = "
            f"{unique_geometric}/{CASE_COUNT}"
        )

        print(
            f"unique exact zero           = "
            f"{unique_zero}/{CASE_COUNT}"
        )

        print()

    print("============================================================")
    print("TEST 2: DETAILED CASE")
    print("============================================================")

    case = cases[0]
    k = 16

    candidates = enumerate_candidates(
        case,
        k,
    )

    geometric = geometric_candidates(
        candidates,
        case["s"],
    )

    zeros = exact_zero(candidates)

    print_case(case, k)

    print()

    print(
        f"all candidates       = "
        f"{len(candidates)}"
    )

    print(
        f"geometric candidates = "
        f"{len(geometric)}"
    )

    print(
        f"exact zero candidates = "
        f"{len(zeros)}"
    )

    print()

    print("Exact integer zeros:")
    if not zeros:
        print("NONE")
    else:
        for candidate in zeros:
            print(
                f"x={candidate['x']} "
                f"y={candidate['y']}"
            )

    print()

    print("First 20 geometric candidates:")

    for candidate in geometric[:20]:
        x = candidate["x"]
        y = candidate["y"]

        p_candidate = (
            case["s"] + 1 - x
        )

        q_candidate = (
            case["s"] + x + y
        )

        print(
            f"x={x:7d} "
            f"y={y:7d} "
            f"p'={p_candidate:7d} "
            f"q'={q_candidate:7d} "
            f"P={candidate['residue']}"
        )

    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT}")


if __name__ == "__main__":
    main()