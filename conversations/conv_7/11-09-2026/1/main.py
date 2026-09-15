import math
import random


EXPERIMENT_NUMBER = 23


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
            is_prime[start : limit + 1 : p] = b"\x00" * count

    return [n for n in range(2, limit + 1) if is_prime[n]]


def polynomial(n, s, d, x, y):
    """
    P_S(x,y) =
        x^2 + (y-1)x + (D-s) - (s+1)y
    """
    return (
        x * x
        + (y - 1) * x
        + (d - s)
        - (s + 1) * y
    )


def binomial_weight(x, y):
    """
    B(x,y) = C(x,y) * 2^(x-y).

    For the binomial transform we only use 0 <= y <= x.
    """
    if y < 0 or y > x:
        return 0

    return math.comb(x, y) * (1 << (x - y))


def direct_binomial_transform(n, s, d, x):
    """
    T(x) = sum_{y=0}^x P(x,y) * C(x,y) * 2^(x-y)
    """
    total = 0

    for y in range(x + 1):
        total += polynomial(n, s, d, x, y) * binomial_weight(x, y)

    return total


def closed_binomial_transform(s, d, x):
    """
    Derived closed form:

    T(x) =
        3^(x-1) *
        [4x^2 - (s+4)x + 3D - 3s]

    """
    if x == 0:
        raise ValueError("This experiment uses x >= 1.")

    h = (
        4 * x * x
        - (s + 4) * x
        + 3 * d
        - 3 * s
    )

    return (3 ** (x - 1)) * h


def transform_core(s, d, x):
    """
    The part of the transform after removing 3^(x-1):

        H(x) = 4x^2 - (s+4)x + 3D - 3s
    """
    return (
        4 * x * x
        - (s + 4) * x
        + 3 * d
        - 3 * s
    )


def construct_coordinates(p, q):
    """Construct the user's factor coordinates."""
    n = p * q
    s = math.isqrt(n)
    d = n - s * s

    x = s + 1 - p
    y = p + q - 2 * s - 1

    return n, s, d, x, y


def test_case(p, q):
    """
    Test one known semiprime.
    """
    n, s, d, x, y = construct_coordinates(p, q)

    p_from_x = s + 1 - x
    q_from_xy = s + x + y

    P_zero = polynomial(n, s, d, x, y)

    direct = direct_binomial_transform(n, s, d, x)
    closed = closed_binomial_transform(s, d, x)

    H = transform_core(s, d, x)

    gcd_h_n = math.gcd(abs(H), n)

    # Algebraic prediction at x = s+1-p:
    #
    # H(x) = p * (4p + 3q - 7s - 4)
    predicted_H = p * (4 * p + 3 * q - 7 * s - 4)

    return {
        "n": n,
        "p": p,
        "q": q,
        "s": s,
        "D": d,
        "x": x,
        "y": y,
        "P": P_zero,
        "recovered_p": p_from_x,
        "recovered_q": q_from_xy,
        "T_direct": direct,
        "T_closed": closed,
        "transform_match": direct == closed,
        "H": H,
        "H_predicted": predicted_H,
        "H_formula_match": H == predicted_H,
        "gcd_H_N": gcd_h_n,
        "gcd_equals_p": gcd_h_n == p,
    }


def generate_semiprimes(primes, count, seed):
    """Generate distinct semiprime factor pairs."""
    random.seed(seed)

    results = set()

    while len(results) < count:
        p = random.choice(primes)
        q = random.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        results.add((p, q))

    return sorted(results)


def run_experiment():
    print(f"START EXPERIMENT {EXPERIMENT_NUMBER}")
    print()
    print("Hypothesis:")
    print("  T(x) = sum_y P(x,y) * C(x,y) * 2^(x-y)")
    print("  T(x) = 3^(x-1) * H(x)")
    print("  H(x) = 4x^2 - (s+4)x + 3D - 3s")
    print()
    print("At the true factor coordinate x = s+1-p:")
    print("  H(x) = p * (4p + 3q - 7s - 4)")
    print()

    primes = sieve_primes(5000)
    cases = generate_semiprimes(primes[:500], 300, seed=23)

    tested = []
    for p, q in cases:
        tested.append(test_case(p, q))

    total = len(tested)

    zero_ok = sum(row["P"] == 0 for row in tested)
    recovery_ok = sum(
        row["recovered_p"] * row["recovered_q"] == row["n"]
        for row in tested
    )
    transform_ok = sum(row["transform_match"] for row in tested)
    formula_ok = sum(row["H_formula_match"] for row in tested)
    gcd_ok = sum(row["gcd_equals_p"] for row in tested)

    print("AGGREGATE RESULTS")
    print("-----------------")
    print(f"Cases tested                    : {total}")
    print(f"P(x,y) = 0                     : {zero_ok}/{total}")
    print(f"Recovered p*q = N              : {recovery_ok}/{total}")
    print(f"Direct transform = closed form : {transform_ok}/{total}")
    print(f"H(x) factor formula            : {formula_ok}/{total}")
    print(f"gcd(H(x), N) = p               : {gcd_ok}/{total}")
    print()

    print("SELECTED CASES")
    print("--------------")
    print(
        "N | p | q | s | D | x | y | "
        "P | H(x) | gcd(H,N)"
    )

    # Pick representative cases from different factor distances.
    selected = []

    for row in tested:
        distance = row["q"] - row["p"]

        if distance < 100:
            selected.append(row)

    selected += [
        row for row in tested
        if 100 <= row["q"] - row["p"] < 1000
    ]

    selected += [
        row for row in tested
        if row["q"] - row["p"] >= 1000
    ]

    # Remove duplicates while preserving order.
    unique = []
    seen = set()

    for row in selected:
        key = row["n"]
        if key not in seen:
            seen.add(key)
            unique.append(row)

    for row in unique[:15]:
        print(
            f"{row['n']} | "
            f"{row['p']} | "
            f"{row['q']} | "
            f"{row['s']} | "
            f"{row['D']} | "
            f"{row['x']} | "
            f"{row['y']} | "
            f"{row['P']} | "
            f"{row['H']} | "
            f"{row['gcd_H_N']}"
        )

    print()
    print("LOCAL CHECKS")
    print("------------")

    # Explicitly check several neighboring x values around the true point.
    for row in tested[:8]:
        n = row["n"]
        p = row["p"]
        q = row["q"]
        s = row["s"]
        d = row["D"]
        x0 = row["x"]

        print(
            f"N={n}, p={p}, q={q}, true_x={x0}"
        )

        for dx in (-2, -1, 0, 1, 2):
            x = x0 + dx

            if x < 1:
                continue

            h = transform_core(s, d, x)
            g = math.gcd(abs(h), n)

            print(
                f"  x={x:4d}  "
                f"H={h:>16d}  "
                f"gcd(H,N)={g}"
            )

        print()

    print("INTERPRETATION")
    print("--------------")
    print("1. The arithmetic-sequence/binomial transform identity is exact.")
    print("2. At the factor coordinate x=s+1-p, H(x) contains p as an exact factor.")
    print("3. gcd(H(x),N) therefore recovers p at the known factor coordinate.")
    print("4. The experiment does NOT yet show how to find x efficiently.")
    print("5. The useful remaining question is whether the binomial/ArithSeq")
    print("   structure gives a way to locate that x without searching factors.")

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT_NUMBER}")


if __name__ == "__main__":
    run_experiment()