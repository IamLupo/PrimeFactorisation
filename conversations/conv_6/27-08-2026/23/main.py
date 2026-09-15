import math
import secrets
import sympy
from collections import Counter, defaultdict


# ============================================================================
# SETTINGS
# ============================================================================

FACTOR_LOW = 10_000
FACTOR_HIGH = 100_000

TRIALS = 300

MODULI = [3, 5, 7, 11, 13, 17, 19, 23]
M = math.prod(MODULI)
STEP = 2 * M

T_MIN = -25
T_MAX = 25


# ============================================================================
# RANDOM PRIME GENERATION
# ============================================================================

def random_prime(low, high):
    while True:
        x = secrets.randbelow(high - low + 1) + low

        if x % 2 == 0:
            x += 1

        if x > high:
            continue

        if sympy.isprime(x):
            return x


# ============================================================================
# UNIQUE SEMIPRIME
# ============================================================================

def generate_unique_semiprime(seen):
    while True:
        p = random_prime(FACTOR_LOW, FACTOR_HIGH)
        q = random_prime(FACTOR_LOW, FACTOR_HIGH)

        p, q = sorted((p, q))
        n = p * q

        if n not in seen:
            seen.add(n)
            return p, q, n


# ============================================================================
# FACTOR PAIRS OF N_t
# ============================================================================

def shifted_prime_pairs(value):
    """
    Return prime factor pairs (a,b) such that

        a*b = value

    and both factors lie inside FACTOR_LOW..FACTOR_HIGH.

    Full factorization is used because value is only around 1e10.
    """

    if value <= 0:
        return []

    factors = sympy.factorint(value)

    # We are specifically interested in semiprimes a*b=value.
    total_factor_count = sum(factors.values())

    if total_factor_count != 2:
        return []

    prime_list = []

    for prime, exponent in factors.items():
        prime_list.extend([prime] * exponent)

    if len(prime_list) != 2:
        return []

    a, b = sorted(prime_list)

    if not (
        FACTOR_LOW <= a <= FACTOR_HIGH
        and FACTOR_LOW <= b <= FACTOR_HIGH
    ):
        return []

    return [(a, b)]


# ============================================================================
# MAIN DATA
# ============================================================================

seen_n = set()

all_events = []
t_distribution = Counter()

trial_count_by_t = Counter()

unique_collision_n = set()

# Store relationships between actual factors and shifted factors.
delta_records = defaultdict(list)

# Store the actual shifted values.
shifted_values = defaultdict(list)


# ============================================================================
# HEADER
# ============================================================================

print("=" * 100)
print("SHIFTED-PRODUCT t-LATTICE EXPERIMENT")
print("=" * 100)

print(f"M           = {M:,}")
print(f"2M          = {STEP:,}")
print(
    f"factor range = "
    f"{FACTOR_LOW:,} - {FACTOR_HIGH:,}"
)
print(
    f"t range      = "
    f"[{T_MIN}, {T_MAX}]"
)
print(f"unique trials = {TRIALS}")
print()


# ============================================================================
# EXPERIMENT
# ============================================================================

for trial in range(1, TRIALS + 1):

    p, q, n = generate_unique_semiprime(seen_n)

    trial_events = []

    for t in range(T_MIN, T_MAX + 1):

        # t=0 is the actual factorization.
        if t == 0:
            continue

        shifted = n + STEP * t

        if shifted <= 0:
            continue

        pairs = shifted_prime_pairs(shifted)

        for a, b in pairs:

            # Sanity check.
            if a * b != shifted:
                raise RuntimeError(
                    "Factorization verification failed"
                )

            dp = a - p
            dq = b - q

            event = {
                "trial": trial,
                "n": n,
                "p": p,
                "q": q,
                "a": a,
                "b": b,
                "t": t,
                "shifted": shifted,
                "dp": dp,
                "dq": dq,
            }

            trial_events.append(event)
            all_events.append(event)

            t_distribution[t] += 1
            trial_count_by_t[t] += 1

            unique_collision_n.add(n)

            delta_records[t].append(
                (dp, dq)
            )

            shifted_values[t].append(
                shifted
            )

    if trial_events:

        print(
            f"trial={trial:3} "
            f"p={p:,} "
            f"q={q:,} "
            f"n={n:,} "
            f"events={len(trial_events)}"
        )

        for event in sorted(
            trial_events,
            key=lambda x: x["t"]
        ):

            print(
                f"    t={event['t']:+3} "
                f"N_t={event['shifted']:,} "
                f"{event['a']:,} * {event['b']:,}"
            )


# ============================================================================
# SUMMARY
# ============================================================================

print()
print("=" * 100)
print("SUMMARY")
print("=" * 100)

print(
    f"unique semiprimes tested = {len(seen_n)}"
)

print(
    f"total shifted semiprime events = "
    f"{len(all_events)}"
)

print(
    f"unique n with events = "
    f"{len(unique_collision_n)}"
)

print()


# ============================================================================
# t DISTRIBUTION
# ============================================================================

print("=" * 100)
print("t DISTRIBUTION")
print("=" * 100)

if not t_distribution:

    print("No shifted semiprime events found.")

else:

    for t in sorted(t_distribution):

        print(
            f"t={t:+4} : "
            f"{t_distribution[t]:4} events "
            f"across "
            f"{trial_count_by_t[t]:4} unique n"
        )


# ============================================================================
# POSITIVE / NEGATIVE COMPARISON
# ============================================================================

print()
print("=" * 100)
print("SYMMETRIC t COMPARISON")
print("=" * 100)

for magnitude in range(1, max(abs(T_MIN), abs(T_MAX)) + 1):

    negative = t_distribution.get(
        -magnitude,
        0
    )

    positive = t_distribution.get(
        +magnitude,
        0
    )

    if negative or positive:

        print(
            f"|t|={magnitude:3} "
            f"negative={negative:4} "
            f"positive={positive:4}"
        )


# ============================================================================
# FACTOR-SHIFT GEOMETRY
# ============================================================================

print()
print("=" * 100)
print("FACTOR-SHIFT GEOMETRY")
print("=" * 100)

for t in sorted(delta_records):

    records = delta_records[t]

    dp_values = [
        x[0]
        for x in records
    ]

    dq_values = [
        x[1]
        for x in records
    ]

    print()
    print(
        f"t={t:+4} "
        f"count={len(records)}"
    )

    print(
        f"    dp range = "
        f"{min(dp_values):+,} ... "
        f"{max(dp_values):+,}"
    )

    print(
        f"    dq range = "
        f"{min(dq_values):+,} ... "
        f"{max(dq_values):+,}"
    )


# ============================================================================
# PRODUCT IDENTITY VERIFICATION
# ============================================================================

print()
print("=" * 100)
print("PRODUCT IDENTITY VERIFICATION")
print("=" * 100)

failures = 0

for event in all_events:

    lhs = (
        event["a"] * event["b"]
        - event["p"] * event["q"]
    )

    rhs = STEP * event["t"]

    if lhs != rhs:

        failures += 1

        print(
            "FAIL:",
            event
        )

print(
    f"verified events = {len(all_events)}"
)

print(
    f"failures        = {failures}"
)


# ============================================================================
# RECONSTRUCT t FROM FACTORS
# ============================================================================

print()
print("=" * 100)
print("t RECONSTRUCTION")
print("=" * 100)

reconstruction_failures = 0

for event in all_events:

    reconstructed = (
        event["a"] * event["b"]
        - event["n"]
    ) // STEP

    if reconstructed != event["t"]:

        reconstruction_failures += 1

        print(
            "FAIL:",
            event
        )

print(
    f"reconstruction failures = "
    f"{reconstruction_failures}"
)


# ============================================================================
# COLLISION PAIR RATIO
# ============================================================================

print()
print("=" * 100)
print("SHIFTED FACTOR NORMALIZATION")
print("=" * 100)

for t in sorted(delta_records):

    records = delta_records[t]

    if not records:
        continue

    # Normalize relative to factor range.
    normalized_dp = [
        dp / FACTOR_HIGH
        for dp, dq in records
    ]

    normalized_dq = [
        dq / FACTOR_HIGH
        for dp, dq in records
    ]

    print()
    print(
        f"t={t:+4}"
    )

    print(
        f"    normalized dp: "
        f"{min(normalized_dp):+.6f} ... "
        f"{max(normalized_dp):+.6f}"
    )

    print(
        f"    normalized dq: "
        f"{min(normalized_dq):+.6f} ... "
        f"{max(normalized_dq):+.6f}"
    )


# ============================================================================
# MOST EXTREME EVENTS
# ============================================================================

print()
print("=" * 100)
print("EXTREME EVENTS")
print("=" * 100)

if all_events:

    largest_abs_t = sorted(
        all_events,
        key=lambda x: abs(x["t"]),
        reverse=True
    )[:20]

    largest_dp = sorted(
        all_events,
        key=lambda x: abs(x["dp"]),
        reverse=True
    )[:20]

    print()
    print("Largest |t|:")

    for event in largest_abs_t:

        print(
            f"  t={event['t']:+4} "
            f"actual=({event['p']:,},{event['q']:,}) "
            f"shifted=({event['a']:,},{event['b']:,})"
        )

    print()
    print("Largest |dp|:")

    for event in largest_dp:

        print(
            f"  t={event['t']:+4} "
            f"dp={event['dp']:+,} "
            f"dq={event['dq']:+,} "
            f"actual=({event['p']:,},{event['q']:,}) "
            f"shifted=({event['a']:,},{event['b']:,})"
        )


# ============================================================================
# UNIQUE SHIFTED VALUES
# ============================================================================

print()
print("=" * 100)
print("UNIQUE SHIFTED SEMIPRIMES")
print("=" * 100)

for t in sorted(shifted_values):

    values = sorted(
        set(shifted_values[t])
    )

    print(
        f"t={t:+4} "
        f"unique N_t values={len(values)}"
    )

    for value in values[:10]:

        print(
            f"    {value:,}"
        )

    if len(values) > 10:

        print(
            f"    ... "
            f"{len(values) - 10} more"
        )


# ============================================================================
# FINAL
# ============================================================================

print()
print("=" * 100)
print("FINAL")
print("=" * 100)

print(
    "Every recorded collision satisfies:"
)

print(
    "    a*b = n + 2*M*t"
)

print(
    "and therefore:"
)

print(
    "    a*b - p*q = 2*M*t"
)

print()
print("Experiment complete.")
