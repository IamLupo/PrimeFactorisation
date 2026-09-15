import math
import random
import time

from sympy import nextprime


EXPERIMENT_ID = 138
CASES = 100


def generate_semiprime():
    while True:
        p = int(nextprime(random.randint(2000, 5000)))

        ratio = random.uniform(1.05, 1.50)
        q = int(nextprime(int(p * ratio)))

        if q > p:
            return p, q, p * q


def initial_interval(n):
    """
    We know:

        p <= sqrt(N)

    and because q/p <= 1.50 in this experiment,

        p >= sqrt(N / 1.50).
    """

    lower = math.ceil(math.sqrt(n / 1.50))
    upper = math.isqrt(n)

    return lower, upper


def make_binary_probe(lower, upper):
    """
    Construct a probe with:

        m = k - 2 = upper
        d = upper - midpoint

    where

        midpoint = floor((lower + upper) / 2).

    Since upper <= floor(sqrt(N)) < q, q cannot divide C(m,d).

    For every p in [lower, upper]:

        floor(m/p) = 1

    and Lucas gives:

        p | C(m,d)  <=>  p > midpoint.
    """

    if lower >= upper:
        return None

    midpoint = (lower + upper) // 2

    m = upper
    d = upper - midpoint
    k = m + 2

    return k, d, midpoint


def evaluate_probe(n, k, d):
    m = k - 2

    if d < 0:
        raise ValueError("d < 0")

    if d >= k - 1:
        return n

    return math.gcd(math.comb(m, d), n)


def classify_signal(g, n):
    if g == 1:
        return "MISS"

    if g == n:
        return "N"

    return "HIT"


def verify_probe_geometry(lower, upper, k, d):
    m = k - 2
    midpoint = m - d

    # d must be strictly smaller than every possible p.
    if d >= lower:
        return False, "d >= lower"

    # m must lie in the quotient-1 region for all possible p.
    if m < upper:
        return False, "m < upper"

    if m >= 2 * lower:
        return False, "m >= 2*lower"

    # m must stay below q.
    # Since q >= sqrt(N), this is checked outside where N is available.
    return True, "OK"


def run_binary_search(n):
    lower, upper = initial_interval(n)

    history = []

    while lower < upper:
        probe = make_binary_probe(lower, upper)

        if probe is None:
            break

        k, d, midpoint = probe

        valid, reason = verify_probe_geometry(
            lower,
            upper,
            k,
            d,
        )

        if not valid:
            return {
                "success": False,
                "reason": reason,
                "history": history,
                "factor_candidate": None,
            }

        g = evaluate_probe(n, k, d)
        signal = classify_signal(g, n)

        history.append(
            {
                "lower": lower,
                "upper": upper,
                "midpoint": midpoint,
                "k": k,
                "d": d,
                "gcd": g,
                "signal": signal,
            }
        )

        if signal == "HIT":
            # p > midpoint
            lower = midpoint + 1

        elif signal == "MISS":
            # p <= midpoint
            upper = midpoint

        else:
            return {
                "success": False,
                "reason": "unexpected N",
                "history": history,
                "factor_candidate": None,
            }

    candidate = lower if lower == upper else None

    return {
        "success": candidate is not None,
        "reason": "OK",
        "history": history,
        "factor_candidate": candidate,
    }


def validate_result(p, result):
    if not result["success"]:
        return False

    return result["factor_candidate"] == p


def print_case(index, p, q, n, result):
    history = result["history"]

    print(f"CASE {index}")
    print(f"p={p}")
    print(f"q={q}")
    print(f"N={n}")

    if history:
        print(
            f"initial_interval="
            f"[{history[0]['lower']},{history[0]['upper']}]"
        )

    print(f"probes={len(history)}")
    print(f"candidate={result['factor_candidate']}")
    print(f"success={result['success']}")
    print(f"reason={result['reason']}")
    print()

    for i, item in enumerate(history, 1):
        print(
            f"  {i:02d}: "
            f"[{item['lower']},{item['upper']}] "
            f"mid={item['midpoint']} "
            f"k={item['k']} "
            f"d={item['d']} "
            f"gcd={item['gcd']} "
            f"{item['signal']}"
        )

    print()


def main():
    print(f"START EXPERIMENT {EXPERIMENT_ID}")
    print()

    random.seed(EXPERIMENT_ID)

    successes = 0
    failures = 0
    total_probes = 0
    max_probes = 0

    start = time.perf_counter()

    for case in range(1, CASES + 1):
        p, q, n = generate_semiprime()

        result = run_binary_search(n)

        if validate_result(p, result):
            successes += 1
        else:
            failures += 1

        probes = len(result["history"])

        total_probes += probes
        max_probes = max(max_probes, probes)

        print_case(
            case,
            p,
            q,
            n,
            result,
        )

    elapsed = time.perf_counter() - start

    print("SUMMARY")
    print(f"cases={CASES}")
    print(f"successes={successes}")
    print(f"failures={failures}")
    print(f"success_rate={successes / CASES:.6f}")
    print(f"average_probes={total_probes / CASES:.3f}")
    print(f"max_probes={max_probes}")
    print(f"runtime_seconds={elapsed:.6f}")
    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT_ID}")


if __name__ == "__main__":
    main()