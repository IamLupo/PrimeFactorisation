import math
import random


EXPERIMENT_NUMBER = 24


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


def polynomial_H(n, s, d, x):
    """
    H(x) = 4x^2 - (s+4)x + 3D - 3s

    where:
        s = floor(sqrt(N))
        D = N - s^2
    """
    return (
        4 * x * x
        - (s + 4) * x
        + 3 * d
        - 3 * s
    )


def delta_H(s, x):
    """
    Forward difference:
        H(x+1) - H(x) = 8x - s
    """
    return 8 * x - s


def expected_factor_coordinate(s, p):
    """Known factor coordinate x = s + 1 - p."""
    return s + 1 - p


def discriminant_H(n, s):
    """
    Discriminant of H:

        Delta_H = (7s+4)^2 - 48N
    """
    return (7 * s + 4) ** 2 - 48 * n


def classify_gcd(g, p, q):
    """Classify a gcd hit."""
    if g == 1:
        return "NONE"
    if g == p and g != q:
        return "P"
    if g == q and g != p:
        return "Q"
    if g == p * q:
        return "N"
    return "OTHER"


def find_gcd_roots(n, p, q, s, d):
    """
    Find every x in [0, s+1] for which gcd(H(x), N) > 1.
    """
    roots = []

    for x in range(0, s + 2):
        h = polynomial_H(n, s, d, x)
        g = math.gcd(abs(h), n)

        if g > 1:
            roots.append(
                {
                    "x": x,
                    "H": h,
                    "gcd": g,
                    "type": classify_gcd(g, p, q),
                }
            )

    return roots


def modular_roots(n, p, q, s, d):
    """
    Find roots of H modulo p and q over the relevant x interval.
    """
    roots_p = []
    roots_q = []

    for x in range(0, s + 2):
        h = polynomial_H(n, s, d, x)

        if h % p == 0:
            roots_p.append(x)

        if h % q == 0:
            roots_q.append(x)

    return roots_p, roots_q


def generate_semiprimes(primes, count, seed):
    """Generate distinct semiprime factor pairs."""
    random.seed(seed)

    pairs = set()

    while len(pairs) < count:
        p = random.choice(primes)
        q = random.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        pairs.add((p, q))

    return sorted(pairs)


def analyze_case(p, q):
    """Analyze one semiprime."""
    n = p * q
    s = math.isqrt(n)
    d = n - s * s

    xp = expected_factor_coordinate(s, p)
    xq = expected_factor_coordinate(s, q)

    roots = find_gcd_roots(n, p, q, s, d)
    roots_p, roots_q = modular_roots(n, p, q, s, d)

    discriminant = discriminant_H(n, s)
    discriminant_square = (
        discriminant >= 0
        and math.isqrt(discriminant) ** 2 == discriminant
    )

    return {
        "n": n,
        "p": p,
        "q": q,
        "s": s,
        "d": d,
        "xp": xp,
        "xq": xq,
        "roots": roots,
        "roots_p": roots_p,
        "roots_q": roots_q,
        "discriminant": discriminant,
        "discriminant_square": discriminant_square,
    }


def run_experiment():
    print(f"START EXPERIMENT {EXPERIMENT_NUMBER}")
    print()

    primes = sieve_primes(5000)

    # Keep the experiment manageable while covering:
    # - close factors
    # - very unequal factors
    # - small factors
    # - larger factors
    cases = generate_semiprimes(primes[:700], 300, seed=24)

    analyses = [
        analyze_case(p, q)
        for p, q in cases
    ]

    print("H(x) GCD-ROOT EXPERIMENT")
    print("------------------------")
    print()
    print("For each semiprime N=p*q:")
    print("  H(x) = 4x^2 - (s+4)x + 3D - 3s")
    print("  Search 0 <= x <= s+1")
    print("  Record every x with gcd(H(x),N) > 1")
    print()

    total = len(analyses)

    cases_with_hits = sum(
        1 for a in analyses
        if a["roots"]
    )

    cases_with_p_hit = sum(
        1 for a in analyses
        if any(r["type"] == "P" for r in a["roots"])
    )

    cases_with_q_hit = sum(
        1 for a in analyses
        if any(r["type"] == "Q" for r in a["roots"])
    )

    cases_with_N_hit = sum(
        1 for a in analyses
        if any(r["type"] == "N" for r in a["roots"])
    )

    xp_is_p_root = sum(
        1 for a in analyses
        if a["xp"] in a["roots_p"]
    )

    xq_is_q_root = sum(
        1 for a in analyses
        if a["xq"] in a["roots_q"]
    )

    discriminant_square_count = sum(
        a["discriminant_square"]
        for a in analyses
    )

    print(f"Cases tested                  : {total}")
    print(f"Cases with gcd root(s)        : {cases_with_hits}/{total}")
    print(f"Cases with P-root hit         : {cases_with_p_hit}/{total}")
    print(f"Cases with Q-root hit         : {cases_with_q_hit}/{total}")
    print(f"Cases with N-root hit         : {cases_with_N_hit}/{total}")
    print(f"x_p is root modulo p          : {xp_is_p_root}/{total}")
    print(f"x_q is root modulo q          : {xq_is_q_root}/{total}")
    print(f"Discriminant is square        : {discriminant_square_count}/{total}")
    print()

    print("CASES WITH MULTIPLE GCD ROOTS")
    print("-----------------------------")

    multiple_cases = [
        a for a in analyses
        if len(a["roots"]) >= 2
    ]

    for a in multiple_cases[:25]:
        root_string = ", ".join(
            f"x={r['x']}:{r['type']}({r['gcd']})"
            for r in a["roots"]
        )

        print(
            f"N={a['n']} "
            f"p={a['p']} "
            f"q={a['q']} "
            f"s={a['s']} "
            f"xp={a['xp']} "
            f"xq={a['xq']} "
            f"roots=[{root_string}]"
        )

    print()
    print("CASES WHERE H(x)=0")
    print("------------------")

    zero_cases = []

    for a in analyses:
        for r in a["roots"]:
            if r["H"] == 0:
                zero_cases.append((a, r))

    for a, r in zero_cases[:25]:
        print(
            f"N={a['n']} "
            f"p={a['p']} "
            f"q={a['q']} "
            f"s={a['s']} "
            f"x={r['x']} "
            f"H=0 "
            f"discriminant={a['discriminant']} "
            f"sqrt_discriminant="
            f"{math.isqrt(a['discriminant']) if a['discriminant'] >= 0 else 'NA'}"
        )

    print()
    print("FIRST ROOT STRUCTURE SAMPLES")
    print("----------------------------")

    shown = 0

    for a in analyses:
        if not a["roots"]:
            continue

        print(
            f"N={a['n']} "
            f"p={a['p']} "
            f"q={a['q']} "
            f"s={a['s']}"
        )

        print(
            f"  expected xp={a['xp']} "
            f"expected xq={a['xq']}"
        )

        print(
            f"  roots modulo p: {a['roots_p']}"
        )

        print(
            f"  roots modulo q: {a['roots_q']}"
        )

        print(
            "  gcd roots: "
            + ", ".join(
                f"x={r['x']} -> gcd={r['gcd']} ({r['type']}), H={r['H']}"
                for r in a["roots"]
            )
        )

        shown += 1

        if shown >= 20:
            break

    print()
    print("DISTANCE FROM KNOWN p-ROOT")
    print("---------------------------")

    distances = []

    for a in analyses:
        xp = a["xp"]

        for r in a["roots"]:
            distances.append(
                (
                    abs(r["x"] - xp),
                    a["n"],
                    a["p"],
                    a["q"],
                    xp,
                    r["x"],
                    r["type"],
                )
            )

    distances.sort()

    for row in distances[:30]:
        distance, n, p, q, xp, x, root_type = row

        print(
            f"N={n} "
            f"p={p} "
            f"q={q} "
            f"xp={xp} "
            f"x={x} "
            f"distance={distance} "
            f"type={root_type}"
        )

    print()
    print("FINITE-DIFFERENCE CHECK")
    print("-----------------------")

    finite_difference_ok = 0

    for a in analyses:
        s = a["s"]
        d = a["d"]

        for x in range(0, min(a["s"] + 1, 20)):
            h0 = polynomial_H(a["n"], s, d, x)
            h1 = polynomial_H(a["n"], s, d, x + 1)

            if h1 - h0 == delta_H(s, x):
                finite_difference_ok += 1

    print(
        "H(x+1)-H(x) = 8x-s checks: "
        f"{finite_difference_ok}"
    )

    print()
    print("INTERPRETATION")
    print("--------------")
    print("1. Map every modular root of H(x), not only x=s+1-p.")
    print("2. Determine whether the q-factor has a systematic x-coordinate.")
    print("3. Determine whether extra gcd roots have a simple relation to xp or xq.")
    print("4. Check whether H(x)=0 cases correspond to a special factor geometry.")
    print("5. The target is a rule predicting a useful root x from N alone.")

    print()
    print(f"FINISHED EXPERIMENT {EXPERIMENT_NUMBER}")


if __name__ == "__main__":
    run_experiment()
