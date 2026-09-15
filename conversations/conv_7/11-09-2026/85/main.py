import math
import random
from sympy import isprime


EXPERIMENT = 106

CASES = 500

P_MIN = 100
P_MAX = 12000
Q_MAX = 16000

MAX_Z = 16384


def gcd_state(g, n, p, q):
    if g == 1:
        return "1"

    if g == n:
        return "n"

    if g == p:
        return "p"

    if g == q:
        return "q"

    return "other"


def value_gcd(n, k, z):
    """
    v_k = -C(k-2, 2z-1) when 2z-1 < k-1.
    Otherwise v_k = 0.
    """
    d = 2 * z - 1

    if d >= k - 1:
        return n

    value = math.comb(k - 2, d) % n
    value = (-value) % n

    return math.gcd(value, n)


def lucas_digits(value, p):
    digits = []

    while value > 0:
        digits.append(value % p)
        value //= p

    if not digits:
        digits.append(0)

    return digits


def lucas_binomial_nonzero(m, d, p):
    """
    Lucas theorem:

        C(m,d) != 0 mod p

    iff every base-p digit of d is <=
    the corresponding digit of m.
    """

    m_digits = lucas_digits(m, p)
    d_digits = lucas_digits(d, p)

    length = max(len(m_digits), len(d_digits))

    for i in range(length):
        md = m_digits[i] if i < len(m_digits) else 0
        dd = d_digits[i] if i < len(d_digits) else 0

        if dd > md:
            return False

    return True


def lucas_predict_factor(p, k, z):
    """
    Predict whether p divides v_k.

    Since

        v_k = -C(k-2, d)

    with d = 2z-1,

    p divides v_k iff Lucas says C(k-2,d)=0 mod p.
    """

    d = 2 * z - 1

    if d >= k - 1:
        return True

    return not lucas_binomial_nonzero(
        k - 2,
        d,
        p,
    )


def powers_between(p, q):
    result = []

    k = 1

    while k <= p:
        k *= 2

    while k < q:
        result.append(k)
        k *= 2

    return result


def generate_semiprime():
    while True:
        p = random.randrange(P_MIN, P_MAX + 1)

        if not isprime(p):
            continue

        q = random.randrange(
            max(p + 1, P_MIN),
            Q_MAX + 1,
        )

        if not isprime(q):
            continue

        if p != q:
            return p, q


def dyadic_points(max_z):
    points = []

    z = 1

    while z <= max_z:
        points.append(z)
        z *= 2

    return points


def run_experiment():
    print("=" * 60)
    print(f"START EXPERIMENT {EXPERIMENT}")
    print("=" * 60)
    print()

    print(f"CASES: {CASES}")
    print(f"P RANGE: {P_MIN}..{P_MAX}")
    print(f"Q MAX: {Q_MAX}")
    print(f"MAX Z: {MAX_Z}")
    print()

    random.seed(106)

    total_k_cases = 0
    total_dyadic_tests = 0

    actual_p_hits = 0
    predicted_p_hits = 0
    prediction_matches = 0
    prediction_mismatches = 0

    dyadic_factor_cases = 0
    dyadic_miss_cases = 0

    clean_regime_tests = 0
    clean_regime_matches = 0

    post_p_tests = 0
    post_p_actual_hits = 0

    mismatch_examples = []
    miss_examples = []

    points = dyadic_points(MAX_Z)

    for case_index in range(1, CASES + 1):
        p, q = generate_semiprime()
        n = p * q

        for k in powers_between(p, q):
            total_k_cases += 1

            for z in points:
                d = 2 * z - 1

                if d >= k - 1:
                    break

                total_dyadic_tests += 1

                actual_g = value_gcd(n, k, z)
                actual_state = gcd_state(
                    actual_g,
                    n,
                    p,
                    q,
                )

                actual_p = (
                    actual_state == "p"
                )

                predicted_p = lucas_predict_factor(
                    p,
                    k,
                    z,
                )

                if actual_p:
                    actual_p_hits += 1

                if predicted_p:
                    predicted_p_hits += 1

                if actual_p == predicted_p:
                    prediction_matches += 1
                else:
                    prediction_mismatches += 1

                    if len(mismatch_examples) < 10:
                        mismatch_examples.append(
                            (
                                case_index,
                                n,
                                p,
                                q,
                                k,
                                z,
                                d,
                                actual_g,
                                actual_state,
                                predicted_p,
                            )
                        )

                if d < p:
                    clean_regime_tests += 1

                    if actual_p == predicted_p:
                        clean_regime_matches += 1

                if d >= p:
                    post_p_tests += 1

                    if actual_p:
                        post_p_actual_hits += 1

            # Did this k have a dyadic p hit at all?
            dyadic_hit = False

            for z in points:
                d = 2 * z - 1

                if d >= k - 1:
                    break

                g = value_gcd(n, k, z)

                if g == p:
                    dyadic_hit = True
                    break

            if dyadic_hit:
                dyadic_factor_cases += 1
            else:
                dyadic_miss_cases += 1

                if len(miss_examples) < 15:
                    miss_examples.append(
                        (
                            case_index,
                            n,
                            p,
                            q,
                            k,
                            (k - 2) % p,
                        )
                    )

    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    print(f"semiprime cases              = {CASES}")
    print(f"(N,k) cases                  = {total_k_cases}")
    print(f"dyadic z tests               = {total_dyadic_tests}")

    print(f"actual p hits                = {actual_p_hits}")
    print(f"Lucas predicted p hits       = {predicted_p_hits}")
    print(f"prediction matches            = {prediction_matches}")
    print(f"prediction mismatches         = {prediction_mismatches}")

    print()
    print(f"dyadic factor cases           = {dyadic_factor_cases}")
    print(f"dyadic miss cases             = {dyadic_miss_cases}")

    print()
    print(f"clean d < p tests             = {clean_regime_tests}")
    print(f"clean matches                 = {clean_regime_matches}")

    print()
    print(f"post-p tests                  = {post_p_tests}")
    print(f"post-p actual p hits          = {post_p_actual_hits}")

    if total_dyadic_tests:
        print(
            f"Lucas accuracy                = "
            f"{100.0 * prediction_matches / total_dyadic_tests:.4f}%"
        )

    if clean_regime_tests:
        print(
            f"clean-regime accuracy         = "
            f"{100.0 * clean_regime_matches / clean_regime_tests:.4f}%"
        )

    print()
    print("-" * 60)
    print("PREDICTION MISMATCHES")
    print("-" * 60)

    if mismatch_examples:
        for (
            case_index,
            n,
            p,
            q,
            k,
            z,
            d,
            g,
            state,
            predicted,
        ) in mismatch_examples:
            print(
                f"CASE={case_index} "
                f"N={n} "
                f"p={p} "
                f"q={q} "
                f"k={k} "
                f"z={z} "
                f"d={d} "
                f"gcd={g} "
                f"state={state} "
                f"lucas_factor={predicted}"
            )
    else:
        print("NONE")

    print()
    print("-" * 60)
    print("DYADIC MISS CASES")
    print("-" * 60)

    for (
        case_index,
        n,
        p,
        q,
        k,
        residue,
    ) in miss_examples:
        print(
            f"CASE={case_index} "
            f"N={n} "
            f"p={p} "
            f"q={q} "
            f"k={k} "
            f"(k-2)%p={residue}"
        )

    print()
    print("=" * 60)
    print(f"FINISHED EXPERIMENT {EXPERIMENT}")
    print("=" * 60)


def main():
    run_experiment()


if __name__ == "__main__":
    main()
