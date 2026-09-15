import sympy
import math
from collections import Counter, defaultdict


# ============================================================
# SETTINGS
# ============================================================

FACTOR_LOW = 10_000
FACTOR_HIGH = 100_000

TRIALS = 100

MODULI = [3, 5, 7, 11, 13, 17, 19, 23]

M = math.prod(MODULI)
STEP = 2 * M

# Only inspect a small t window first.
T_MIN = -100
T_MAX = 100


# ============================================================
# GENERATE SEMIPRIME
# ============================================================

def generate_semiprime():
    p = sympy.randprime(
        FACTOR_LOW,
        FACTOR_HIGH
    )

    q = sympy.randprime(
        FACTOR_LOW,
        FACTOR_HIGH
    )

    return p, q, p * q


# ============================================================
# CHECK WHETHER A NUMBER HAS A PRIME PAIR IN RANGE
# ============================================================

def check_factorization(value):
    """
    Factor value and return a pair only if value is exactly
    a product of two primes in the requested range.
    """

    factors = sympy.factorint(value)

    # Need exactly two prime factors counting multiplicity.
    expanded = []

    for prime, exponent in factors.items():
        expanded.extend([prime] * exponent)

    if len(expanded) != 2:
        return None

    a, b = sorted(expanded)

    if not (
        FACTOR_LOW <= a <= FACTOR_HIGH
        and FACTOR_LOW <= b <= FACTOR_HIGH
    ):
        return None

    return a, b


# ============================================================
# MAIN
# ============================================================

print("=" * 100)
print("SHIFTED-PRODUCT / t-LATTICE EXPERIMENT")
print("=" * 100)

print(f"M    = {M:,}")
print(f"2M   = {STEP:,}")
print()

print(
    f"Factor range = "
    f"{FACTOR_LOW:,} - {FACTOR_HIGH:,}"
)

print(
    f"t range = "
    f"[{T_MIN}, {T_MAX}]"
)

print(
    f"Trials = {TRIALS}"
)

print()


all_events = []

t_counter = Counter()


for trial in range(1, TRIALS + 1):

    p, q, n = generate_semiprime()

    actual = tuple(sorted((p, q)))

    events = []

    # --------------------------------------------------------
    # Scan the arithmetic progression
    #
    # N_t = n + 2*M*t
    # --------------------------------------------------------

    for t in range(T_MIN, T_MAX + 1):

        if t == 0:
            continue

        value = n + STEP * t

        if value <= 0:
            continue

        pair = check_factorization(value)

        if pair is None:
            continue

        # Ignore the original factorization if it somehow
        # occurs at t=0 (normally excluded above anyway).
        if pair == actual:
            continue

        a, b = pair

        # Verify the defining identity.
        assert a * b == value
        assert a * b - n == STEP * t

        event = {
            "trial": trial,
            "p": p,
            "q": q,
            "n": n,
            "t": t,
            "a": a,
            "b": b,
            "delta_p": a - actual[0],
            "delta_q": b - actual[1],
        }

        events.append(event)
        all_events.append(event)
        t_counter[t] += 1

    if events:

        print(
            f"trial {trial:3}: "
            f"p={p:6,} "
            f"q={q:6,} "
            f"collisions={len(events)}"
        )

        for e in events:

            print(
                f"    t={e['t']:+4} "
                f"{e['a']:,} * {e['b']:,}"
            )

    elif trial % 10 == 0:

        print(
            f"trial {trial:3}: no shifted "
            f"semiprime collisions"
        )


# ============================================================
# SUMMARY
# ============================================================

print()
print("=" * 100)
print("SUMMARY")
print("=" * 100)

print(
    f"Total shifted semiprime events = "
    f"{len(all_events)}"
)

print()

if all_events:

    print(
        f"Unique trials with events = "
        f"{len(set(e['trial'] for e in all_events))}"
    )

    print()

    print("t distribution:")

    for t, count in sorted(t_counter.items()):

        print(
            f"  t={t:+4} : {count:4}"
        )

else:

    print("No events found.")


# ============================================================
# SYMMETRIC t DISTRIBUTION
# ============================================================

print()
print("=" * 100)
print("SYMMETRIC t COMPARISON")
print("=" * 100)

for t in range(1, max(
    abs(x)
    for x in t_counter
) + 1):

    negative = t_counter.get(-t, 0)
    positive = t_counter.get(+t, 0)

    if negative or positive:

        print(
            f"|t|={t:3}: "
            f"negative={negative:4} "
            f"positive={positive:4}"
        )


# ============================================================
# SMALL |t|
# ============================================================

print()
print("=" * 100)
print("SMALLEST |t| EVENTS")
print("=" * 100)

for e in sorted(
    all_events,
    key=lambda x: (abs(x["t"]), x["trial"])
)[:50]:

    print(
        f"trial={e['trial']:3} "
        f"t={e['t']:+3} "
        f"n={e['n']:,}"
    )

    print(
        f"    actual = "
        f"{e['p']:,} * {e['q']:,}"
    )

    print(
        f"    shifted = "
        f"{e['a']:,} * {e['b']:,}"
    )

    print(
        f"    dp={e['delta_p']:+,} "
        f"dq={e['delta_q']:+,}"
    )


# ============================================================
# CHECK FOR REPEATED FACTOR TRANSFORMATIONS
# ============================================================

print()
print("=" * 100)
print("REPEATED DELTA PATTERNS")
print("=" * 100)

delta_counter = Counter(
    (
        e["delta_p"],
        e["delta_q"]
    )
    for e in all_events
)

for (dp, dq), count in delta_counter.most_common(20):

    if count > 1:

        print(
            f"dp={dp:+,} "
            f"dq={dq:+,} "
            f"count={count}"
        )


# ============================================================
# CHECK FOR LINEAR RELATIONSHIPS
# ============================================================

print()
print("=" * 100)
print("DELTA RELATIONSHIPS")
print("=" * 100)

for e in sorted(
    all_events,
    key=lambda x: abs(x["t"])
)[:100]:

    dp = e["delta_p"]
    dq = e["delta_q"]

    if dp == 0:
        ratio = "undefined"
    else:
        ratio = f"{dq / dp:.8f}"

    print(
        f"t={e['t']:+3} "
        f"dp={dp:+8,} "
        f"dq={dq:+8,} "
        f"dq/dp={ratio}"
    )
