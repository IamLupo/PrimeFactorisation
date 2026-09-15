import math
import random


EXPERIMENT_NUMBER = 36


def is_probable_prime(n):
    """Deterministic Miller-Rabin for 64-bit integers."""
    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small_primes:
        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    bases = (
        2,
        325,
        9375,
        28178,
        450775,
        9780504,
        1795265022,
    )

    for a in bases:
        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        composite = True

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                composite = False
                break

        if composite:
            return False

    return True


def random_prime(bits, rng):
    """Generate a prime with exactly 'bits' bits."""
    while True:
        candidate = rng.getrandbits(bits)

        candidate |= 1 << (bits - 1)
        candidate |= 1

        if is_probable_prime(candidate):
            return candidate


def second_branch_root(s, q):
    """Solve 4*x2 + 3*s = 0 (mod q)."""
    return (-3 * s * pow(4, q - 2, q)) % q


def branch_multiplier(s, q, x2):
    """Compute t from 4*x2 + 3*s = t*q."""
    value = 4 * x2 + 3 * s

    if value % q != 0:
        raise ValueError("Invalid second-branch root.")

    return value // q


def construct_case(p, q):
    """Construct all factor-coordinate quantities."""
    n = p * q
    s = math.isqrt(n)

    x2 = second_branch_root(s, q)
    t = branch_multiplier(s, q, x2)

    xp = s + 1 - p
    y = p + q - 2 * s - 1

    return {
        "N": n,
        "p": p,
        "q": q,
        "s": s,
        "x2": x2,
        "t": t,
        "xp": xp,
        "y": y,
    }


def generate_cases(bits_p, bits_q, count, rng):
    """Generate distinct semiprimes."""
    cases = []
    seen = set()

    while len(cases) < count:
        p = random_prime(bits_p, rng)
        q = random_prime(bits_q, rng)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n in seen:
            continue

        seen.add(n)
        cases.append(construct_case(p, q))

    return cases


def deterministic_groups(cases, key_function, target_ts=(3, 5)):
    """
    Group cases by a key and determine whether the key uniquely
    identifies t among the selected target t values.
    """
    groups = {}

    for case in cases:
        if case["t"] not in target_ts:
            continue

        key = key_function(case)

        groups.setdefault(key, set()).add(case["t"])

    deterministic = 0
    ambiguous = 0

    for values in groups.values():
        if len(values) == 1:
            deterministic += 1
        else:
            ambiguous += 1

    return groups, deterministic, ambiguous


def collision_pairs(cases, key_function, target_ts=(3, 5)):
    """
    Count 3<->5 collisions under a key.
    """
    groups = {}

    for case in cases:
        if case["t"] not in target_ts:
            continue

        key = key_function(case)
        groups.setdefault(key, set()).add(case["t"])

    collisions = []

    for key, values in groups.items():
        if len(values) > 1:
            collisions.append((key, sorted(values)))

    return collisions


def run_experiment():
    print(f"START EXPERIMENT {EXPERIMENT_NUMBER}")
    print()

    print("x_p / y MODULAR BRANCH-DISAMBIGUATION EXPERIMENT")
    print("------------------------------------------------")
    print()
    print("Known relations:")
    print("  x_p = s + 1 - p")
    print("  y   = p + q - 2s - 1")
    print()
    print("For the q-branch:")
    print("  t=3 -> q =  s (mod 4)")
    print("  t=5 -> q = -s (mod 4)")
    print()
    print("Question:")
    print("  Can x_p, y, s, or simple combinations recover")
    print("  the missing q mod 4 information?")
    print()

    rng = random.Random(36)

    scales = [
        (18, 18, 3000),
        (22, 22, 3000),
        (26, 26, 3000),
        (30, 30, 3000),
    ]

    all_cases = []

    for bits_p, bits_q, count in scales:
        cases = generate_cases(
            bits_p,
            bits_q,
            count,
            rng,
        )

        all_cases.extend(cases)

    total = len(all_cases)

    target_cases = [
        case
        for case in all_cases
        if case["t"] in (3, 5)
    ]

    print("AGGREGATE RESULTS")
    print("-----------------")
    print(
        f"Total cases: {total}"
    )
    print(
        f"t=3 or t=5 cases: "
        f"{len(target_cases)}"
    )
    print()

    # ---------------------------------------------------------
    # Verify the basic q mod 4 identities.
    # ---------------------------------------------------------

    print("q MOD 4 BRANCH CHECK")
    print("--------------------")

    t3_ok = 0
    t5_ok = 0

    for case in target_cases:
        if case["t"] == 3:
            if case["q"] % 4 == case["s"] % 4:
                t3_ok += 1

        elif case["t"] == 5:
            if case["q"] % 4 == (-case["s"]) % 4:
                t5_ok += 1

    t3_total = sum(
        case["t"] == 3
        for case in target_cases
    )

    t5_total = sum(
        case["t"] == 5
        for case in target_cases
    )

    print(
        f"t=3: q=s mod4      {t3_ok}/{t3_total}"
    )

    print(
        f"t=5: q=-s mod4     {t5_ok}/{t5_total}"
    )

    print()

    # ---------------------------------------------------------
    # Test simple fingerprints.
    # ---------------------------------------------------------

    fingerprints = {
        "(s mod4)": lambda a: (
            a["s"] % 4,
        ),

        "(x_p mod4)": lambda a: (
            a["xp"] % 4,
        ),

        "(y mod4)": lambda a: (
            a["y"] % 4,
        ),

        "(s,x_p) mod4": lambda a: (
            a["s"] % 4,
            a["xp"] % 4,
        ),

        "(s,y) mod4": lambda a: (
            a["s"] % 4,
            a["y"] % 4,
        ),

        "(x_p,y) mod4": lambda a: (
            a["xp"] % 4,
            a["y"] % 4,
        ),

        "(s,x_p,y) mod4": lambda a: (
            a["s"] % 4,
            a["xp"] % 4,
            a["y"] % 4,
        ),

        "(x_p+y) mod4": lambda a: (
            (a["xp"] + a["y"]) % 4,
        ),

        "(x_p-3y) mod4": lambda a: (
            (a["xp"] - 3 * a["y"]) % 4,
        ),

        "(y-s) mod4": lambda a: (
            (a["y"] - a["s"]) % 4,
        ),
    }

    print("SIMPLE MOD-4 FINGERPRINTS")
    print("------------------------")

    for name, key_function in fingerprints.items():
        groups, deterministic, ambiguous = deterministic_groups(
            target_cases,
            key_function,
        )

        print(
            f"{name:22s} "
            f"groups={len(groups):3d} "
            f"det={deterministic:3d} "
            f"ambiguous={ambiguous:3d}"
        )

    print()

    # ---------------------------------------------------------
    # Direct test:
    #
    # Can the fingerprint distinguish t for each fixed s mod 4?
    # ---------------------------------------------------------

    print("t=3 vs t=5 BY (s mod4, x_p mod4)")
    print("----------------------------------")

    groups = {}

    for case in target_cases:
        key = (
            case["s"] % 4,
            case["xp"] % 4,
        )

        groups.setdefault(key, set()).add(
            case["t"]
        )

    for key in sorted(groups):
        print(
            f"smod4={key[0]} "
            f"x_pmod4={key[1]} "
            f"possible_t={sorted(groups[key])}"
        )

    print()

    # ---------------------------------------------------------
    # Test whether x_p + y reveals p+q relation modulo 4.
    # ---------------------------------------------------------

    print("COMBINATION IDENTITIES")
    print("----------------------")

    combinations = {
        "xp+y": lambda a: (
            a["xp"] + a["y"]
        ),

        "xp-y": lambda a: (
            a["xp"] - a["y"]
        ),

        "xp+3y": lambda a: (
            a["xp"] + 3 * a["y"]
        ),

        "xp-3y": lambda a: (
            a["xp"] - 3 * a["y"]
        ),

        "s+xp+y": lambda a: (
            a["s"] + a["xp"] + a["y"]
        ),

        "s-xp+y": lambda a: (
            a["s"] - a["xp"] + a["y"]
        ),

        "2s+xp+y": lambda a: (
            2 * a["s"] + a["xp"] + a["y"]
        ),
    }

    for name, expression in combinations.items():
        for modulus in [4, 8, 16]:
            groups = {}

            for case in target_cases:
                key = expression(case) % modulus
                groups.setdefault(key, set()).add(
                    case["t"]
                )

            ambiguous = sum(
                len(values) > 1
                for values in groups.values()
            )

            print(
                f"{name:12s} mod{modulus:2d}: "
                f"groups={len(groups):3d} "
                f"ambiguous={ambiguous:3d}"
            )

        print()

    # ---------------------------------------------------------
    # Compare the coordinate fingerprint against q mod 4.
    # ---------------------------------------------------------

    print("CAN (x_p,y,s) RECONSTRUCT q MOD 4?")
    print("----------------------------------")

    for modulus in [4, 8, 16, 32, 64]:
        groups = {}

        for case in target_cases:
            key = (
                case["s"] % modulus,
                case["xp"] % modulus,
                case["y"] % modulus,
            )

            groups.setdefault(key, set()).add(
                case["q"] % 4
            )

        deterministic = sum(
            len(values) == 1
            for values in groups.values()
        )

        print(
            f"mod={modulus:2d}: "
            f"deterministic qmod4 groups="
            f"{deterministic}/{len(groups)}"
        )

    print()

    # ---------------------------------------------------------
    # Find explicit collisions at mod 4.
    # ---------------------------------------------------------

    print("EXPLICIT (x_p,y,s) MOD-4 COLLISIONS")
    print("-----------------------------------")

    key_function = lambda a: (
        a["s"] % 4,
        a["xp"] % 4,
        a["y"] % 4,
    )

    collisions = collision_pairs(
        target_cases,
        key_function,
    )

    print(
        f"Conflicting fingerprints: "
        f"{len(collisions)}"
    )

    shown = 0

    for key, values in collisions:
        print(
            f"fingerprint={key} "
            f"possible_t={values}"
        )

        matching = [
            case
            for case in target_cases
            if key_function(case) == key
        ]

        for case in matching[:6]:
            print(
                f"  N={case['N']} "
                f"p={case['p']} "
                f"q={case['q']} "
                f"s={case['s']} "
                f"x_p={case['xp']} "
                f"y={case['y']} "
                f"t={case['t']}"
            )

        print()

        shown += 1

        if shown >= 20:
            break

    print()

    # ---------------------------------------------------------
    # Look for a direct formula for q mod 4.
    #
    # Test all affine combinations
    #
    #   a*s + b*xp + c*y + d (mod 4)
    #
    # with coefficients 0..3.
    # ---------------------------------------------------------

    print("AFFINE MOD-4 SEARCH")
    print("-------------------")
    print(
        "Searching:"
    )
    print(
        "  f = a*s + b*x_p + c*y + d (mod 4)"
    )
    print(
        "such that f uniquely predicts q mod 4."
    )
    print()

    successful = []

    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    mapping = {}
                    valid = True

                    for case in target_cases:
                        value = (
                            a * case["s"]
                            + b * case["xp"]
                            + c * case["y"]
                            + d
                        ) % 4

                        q_mod4 = case["q"] % 4

                        previous = mapping.get(value)

                        if previous is None:
                            mapping[value] = q_mod4

                        elif previous != q_mod4:
                            valid = False
                            break

                    if valid:
                        successful.append(
                            (a, b, c, d)
                        )

    print(
        f"Successful affine predictors: "
        f"{len(successful)}"
    )

    for coefficients in successful[:30]:
        print(
            f"a={coefficients[0]} "
            f"b={coefficients[1]} "
            f"c={coefficients[2]} "
            f"d={coefficients[3]}"
        )

    print()

    # ---------------------------------------------------------
    # Direct t predictor using affine expressions.
    # ---------------------------------------------------------

    print("AFFINE MOD-4 SEARCH FOR t CLASS")
    print("--------------------------------")

    successful_t = []

    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    mapping = {}
                    valid = True

                    for case in target_cases:
                        value = (
                            a * case["s"]
                            + b * case["xp"]
                            + c * case["y"]
                            + d
                        ) % 4

                        previous = mapping.get(value)

                        if previous is None:
                            mapping[value] = case["t"]

                        elif previous != case["t"]:
                            valid = False
                            break

                    if valid:
                        successful_t.append(
                            (a, b, c, d)
                        )

    print(
        f"Successful affine t predictors: "
        f"{len(successful_t)}"
    )

    for coefficients in successful_t[:30]:
        print(
            f"a={coefficients[0]} "
            f"b={coefficients[1]} "
            f"c={coefficients[2]} "
            f"d={coefficients[3]}"
        )

    print()

    print("INTERPRETATION")
    print("--------------")
    print("1. First check whether (s,x_p,y) mod 4 separates t=3 and t=5.")
    print("2. Then check whether these coordinates reconstruct q mod 4.")
    print("3. The affine search tests whether a simple linear invariant")
    print("   of the arithmetic-sequence coordinates contains the missing")
    print("   branch information.")
    print("4. A successful affine predictor would be substantially more")
    print("   interesting than the raw 2-adic fingerprint.")
    print()

    print(f"FINISHED EXPERIMENT {EXPERIMENT_NUMBER}")


if __name__ == "__main__":
    run_experiment()
