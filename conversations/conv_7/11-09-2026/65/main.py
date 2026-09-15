import math
import random


EXPERIMENT = 83


def is_prime(n):
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False

    r = int(math.isqrt(n))
    d = 3

    while d <= r:
        if n % d == 0:
            return False
        d += 2

    return True


def generate_semiprime(min_prime=1000, max_prime=50000):
    while True:
        p = random.randrange(min_prime, max_prime)
        q = random.randrange(min_prime, max_prime)

        if p > q:
            p, q = q, p

        if p == q:
            continue

        if not is_prime(p) or not is_prime(q):
            continue

        return p, q


def construct_parameters(p, q):
    n = p * q
    s = math.isqrt(n)
    d = n - s * s

    x = s + 1 - p
    y = p + q - 2 * s - 1

    return n, s, d, x, y


def mod_inverse(a, m):
    return pow(a % m, -1, m)


def orbit_point(x, s, d, m):
    denominator = (x - s - 1) % m

    if math.gcd(denominator, m) != 1:
        return None

    numerator = (-((x * x) - x + d - s)) % m

    y = (numerator * mod_inverse(denominator, m)) % m
    x2 = (1 - y - x) % m

    return x2


def is_admissible(x, s, m):
    denominator = (x - s - 1) % m
    return math.gcd(denominator, m) == 1


def orbit_signature(x, s, d, m):
    if not is_admissible(x, s, m):
        return None

    x2 = orbit_point(x, s, d, m)

    if x2 is None:
        return None

    fixed = (x2 == (x % m))

    return {
        "start": x % m,
        "partner": x2,
        "fixed": fixed,
        "orbit_size": 1 if fixed else 2,
    }


def build_admissible_signatures(s, d, m):
    signatures = []

    for x in range(m):
        sig = orbit_signature(x, s, d, m)

        if sig is not None:
            signatures.append(sig)

    return signatures


def true_signature(x_true, s, d, moduli):
    result = []

    for m in moduli:
        sig = orbit_signature(x_true, s, d, m)

        if sig is None:
            result.append((m, None))
        else:
            result.append(
                (
                    m,
                    (
                        sig["fixed"],
                        sig["orbit_size"],
                    ),
                )
            )

    return tuple(result)


def random_signature(s, d, moduli):
    result = []

    for m in moduli:
        candidates = [
            x
            for x in range(m)
            if is_admissible(x, s, m)
        ]

        if not candidates:
            result.append((m, None))
            continue

        x = random.choice(candidates)
        sig = orbit_signature(x, s, d, m)

        if sig is None:
            result.append((m, None))
        else:
            result.append(
                (
                    m,
                    (
                        sig["fixed"],
                        sig["orbit_size"],
                    ),
                )
            )

    return tuple(result)


def centered_candidates(s, window):
    start = max(0, s - window)
    end = s + 1 + window

    return list(range(start, end))


def residue_signature(x, s, d, moduli):
    result = []

    for m in moduli:
        xm = x % m
        sig = orbit_signature(xm, s, d, m)

        if sig is None:
            result.append((m, None))
        else:
            result.append(
                (
                    m,
                    (
                        sig["fixed"],
                        sig["orbit_size"],
                    ),
                )
            )

    return tuple(result)


def count_signature_collisions(x_true, s, d, moduli, window):
    target = residue_signature(x_true, s, d, moduli)

    candidates = centered_candidates(s, window)

    matches = 0
    admissible = 0

    for x in candidates:
        sig = residue_signature(x, s, d, moduli)

        if all(entry[1] is not None for entry in sig):
            admissible += 1

        if sig == target:
            matches += 1

    return matches, admissible, len(candidates)


def run_experiment(num_cases=500):
    moduli = [8, 3, 5, 7, 11, 13, 17, 19, 23]
    window = 200

    exact_signature_matches = 0
    fewer_than_ten_matches = 0

    total_centered_matches = 0
    total_centered_candidates = 0
    total_admissible_candidates = 0

    fixed_true_count = 0
    nonfixed_true_count = 0

    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print("MODULI:")
    print(moduli)
    print()

    print(f"CASES: {num_cases}")
    print(f"CENTERED WINDOW: +/-{window}")
    print()

    for case in range(1, num_cases + 1):
        p, q = generate_semiprime()
        n, s, d, x_true, y_true = construct_parameters(p, q)

        target = true_signature(x_true, s, d, moduli)

        random_matches = 0
        random_trials = 100

        for _ in range(random_trials):
            candidate = random.randrange(0, n)

            if residue_signature(candidate, s, d, moduli) == target:
                random_matches += 1

        matches, admissible, candidate_count = count_signature_collisions(
            x_true,
            s,
            d,
            moduli,
            window,
        )

        total_centered_matches += matches
        total_centered_candidates += candidate_count
        total_admissible_candidates += admissible

        if random_matches > 0:
            exact_signature_matches += 1

        if matches <= 10:
            fewer_than_ten_matches += 1

        true_fixed_count = 0

        for m in moduli:
            sig = orbit_signature(x_true, s, d, m)

            if sig is not None and sig["fixed"]:
                true_fixed_count += 1

        if true_fixed_count > len(moduli) // 2:
            fixed_true_count += 1
        else:
            nonfixed_true_count += 1

        if case <= 10:
            print(
                f"CASE {case}: "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"s={s} "
                f"x_p={x_true} "
                f"y={y_true}"
            )

            print(
                f"  CENTERED MATCHES={matches} "
                f"ADMISSIBLE={admissible}/{candidate_count}"
            )

            print(
                f"  RANDOM SIGNATURE MATCHES="
                f"{random_matches}/{random_trials}"
            )

            print(
                f"  TRUE SIGNATURE={target}"
            )

            print()

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)

    centered_density = (
        total_centered_matches / total_centered_candidates
        if total_centered_candidates
        else 0.0
    )

    admissible_density = (
        total_admissible_candidates / total_centered_candidates
        if total_centered_candidates
        else 0.0
    )

    print(f"CASES: {num_cases}")
    print(
        "CASES WITH RANDOM SIGNATURE COLLISION: "
        f"{exact_signature_matches}/{num_cases}"
    )
    print(
        "CASES WITH <=10 CENTERED SIGNATURE MATCHES: "
        f"{fewer_than_ten_matches}/{num_cases}"
    )
    print(
        "TOTAL CENTERED MATCHES: "
        f"{total_centered_matches}"
    )
    print(
        "TOTAL CENTERED CANDIDATES: "
        f"{total_centered_candidates}"
    )
    print(
        "TOTAL ADMISSIBLE CANDIDATES: "
        f"{total_admissible_candidates}"
    )
    print(
        "CENTERED MATCH DENSITY: "
        f"{centered_density:.8f}"
    )
    print(
        "ADMISSIBLE DENSITY: "
        f"{admissible_density:.8f}"
    )
    print(
        "TRUE ORBIT MORE-OFTEN-FIXED CASES: "
        f"{fixed_true_count}"
    )
    print(
        "TRUE ORBIT MORE-OFTEN-NONFIXED CASES: "
        f"{nonfixed_true_count}"
    )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    random.seed(830011)

    run_experiment(
        num_cases=500
    )


if __name__ == "__main__":
    main()