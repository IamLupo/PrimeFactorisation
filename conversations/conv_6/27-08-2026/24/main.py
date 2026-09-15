import math
import secrets
import sympy
from collections import Counter


# =============================================================================
# SETTINGS
# =============================================================================

FACTOR_LOW = 10_000
FACTOR_HIGH = 100_000

TRIALS = 300

T_MIN = -25
T_MAX = 25

MODULI = [3, 5, 7, 11, 13, 17, 19, 23]

M = math.prod(MODULI)
STEP = 2 * M


# =============================================================================
# PRIME GENERATION
# =============================================================================

def random_prime(low, high):
    while True:
        x = secrets.randbelow(high - low + 1) + low

        if x % 2 == 0:
            x += 1

        if x > high:
            continue

        if sympy.isprime(x):
            return x


# =============================================================================
# ACTUAL SEMIPRIME
# =============================================================================

def generate_unique_semiprime(seen):
    while True:
        p = random_prime(FACTOR_LOW, FACTOR_HIGH)
        q = random_prime(FACTOR_LOW, FACTOR_HIGH)

        p, q = sorted((p, q))
        n = p * q

        if n not in seen:
            seen.add(n)
            return p, q, n


# =============================================================================
# RANDOM CONTROL
# =============================================================================

def generate_control_n(low_n, high_n):
    """
    Generate a random odd integer n0 satisfying

        gcd(n0, 2M) = 1

    so the control has the same basic small-prime exclusion as pq.
    """

    while True:
        n = secrets.randbelow(high_n - low_n + 1) + low_n

        if n % 2 == 0:
            continue

        if math.gcd(n, STEP) != 1:
            continue

        return n


# =============================================================================
# TARGET-RANGE SEMIPRIME TEST
# =============================================================================

def target_semiprime_pair(value):
    """
    Return (a,b) when

        value = a*b

    with both a,b prime and inside the target range.

    Otherwise return None.
    """

    if value <= 0:
        return None

    factors = sympy.factorint(value)

    # Exactly two prime factors counted with multiplicity.
    if sum(factors.values()) != 2:
        return None

    primes = []

    for prime, exponent in factors.items():
        primes.extend([prime] * exponent)

    if len(primes) != 2:
        return None

    a, b = sorted(primes)

    if not (
        FACTOR_LOW <= a <= FACTOR_HIGH
        and FACTOR_LOW <= b <= FACTOR_HIGH
    ):
        return None

    return a, b


# =============================================================================
# STORAGE
# =============================================================================

actual_events = []
control_events = []

actual_t = Counter()
control_t = Counter()

actual_trial_counts = Counter()
control_trial_counts = Counter()

actual_n_with_event = set()
control_n_with_event = set()

seen_actual = set()


# =============================================================================
# HEADER
# =============================================================================

print("=" * 100)
print("ACTUAL vs RANDOM CONTROL t-LATTICE EXPERIMENT")
print("=" * 100)

print(f"M              = {M:,}")
print(f"2M             = {STEP:,}")
print(
    f"factor range   = "
    f"{FACTOR_LOW:,} - {FACTOR_HIGH:,}"
)
print(
    f"t range        = "
    f"[{T_MIN}, {T_MAX}]"
)
print(f"trials          = {TRIALS}")
print()

print(
    "CONTROL:"
)
print(
    "random odd n with gcd(n, 2M) = 1"
)
print()


# =============================================================================
# GENERATE SAME-SCALE CONTROL RANGE
# =============================================================================

MIN_N = FACTOR_LOW * FACTOR_LOW
MAX_N = FACTOR_HIGH * FACTOR_HIGH


# =============================================================================
# MAIN EXPERIMENT
# =============================================================================

for trial in range(1, TRIALS + 1):

    # -------------------------------------------------------------------------
    # REAL SEMIPRIME
    # -------------------------------------------------------------------------

    p, q, n = generate_unique_semiprime(seen_actual)

    trial_actual = []

    for t in range(T_MIN, T_MAX + 1):

        if t == 0:
            continue

        shifted = n + STEP * t

        pair = target_semiprime_pair(shifted)

        if pair is None:
            continue

        a, b = pair

        event = {
            "trial": trial,
            "n": n,
            "p": p,
            "q": q,
            "t": t,
            "shifted": shifted,
            "a": a,
            "b": b,
        }

        actual_events.append(event)
        trial_actual.append(event)

        actual_t[t] += 1
        actual_trial_counts[t] += 1
        actual_n_with_event.add(n)

    # -------------------------------------------------------------------------
    # CONTROL INTEGER
    # -------------------------------------------------------------------------

    control_n = generate_control_n(
        MIN_N,
        MAX_N
    )

    trial_control = []

    for t in range(T_MIN, T_MAX + 1):

        if t == 0:
            continue

        shifted = control_n + STEP * t

        pair = target_semiprime_pair(shifted)

        if pair is None:
            continue

        a, b = pair

        event = {
            "trial": trial,
            "n": control_n,
            "t": t,
            "shifted": shifted,
            "a": a,
            "b": b,
        }

        control_events.append(event)
        trial_control.append(event)

        control_t[t] += 1
        control_trial_counts[t] += 1
        control_n_with_event.add(control_n)

    # -------------------------------------------------------------------------
    # PRINT EVENTS
    # -------------------------------------------------------------------------

    if trial_actual or trial_control:

        print(
            f"trial={trial:3} "
            f"actual_events={len(trial_actual):2} "
            f"control_events={len(trial_control):2}"
        )

        for e in sorted(
            trial_actual,
            key=lambda x: x["t"]
        ):
            print(
                f"    ACTUAL  "
                f"t={e['t']:+3} "
                f"N_t={e['shifted']:,} "
                f"{e['a']:,} * {e['b']:,}"
            )

        for e in sorted(
            trial_control,
            key=lambda x: x["t"]
        ):
            print(
                f"    CONTROL "
                f"t={e['t']:+3} "
                f"N_t={e['shifted']:,} "
                f"{e['a']:,} * {e['b']:,}"
            )


# =============================================================================
# BASIC SUMMARY
# =============================================================================

print()
print("=" * 100)
print("SUMMARY")
print("=" * 100)

print(
    f"actual semiprimes tested       = {TRIALS}"
)

print(
    f"actual collision events        = "
    f"{len(actual_events)}"
)

print(
    f"actual n with >=1 event        = "
    f"{len(actual_n_with_event)}"
)

print()

print(
    f"control integers tested        = {TRIALS}"
)

print(
    f"control collision events       = "
    f"{len(control_events)}"
)

print(
    f"control n with >=1 event       = "
    f"{len(control_n_with_event)}"
)


# =============================================================================
# EVENT RATE
# =============================================================================

actual_rate = len(actual_events) / TRIALS
control_rate = len(control_events) / TRIALS

print()
print("=" * 100)
print("EVENT RATE")
print("=" * 100)

print(
    f"actual events / trial  = "
    f"{actual_rate:.6f}"
)

print(
    f"control events / trial = "
    f"{control_rate:.6f}"
)

if control_rate > 0:

    ratio = actual_rate / control_rate

    print(
        f"actual/control ratio  = "
        f"{ratio:.6f}"
    )

else:

    print(
        "actual/control ratio  = undefined "
        "(zero control events)"
    )


# =============================================================================
# UNIQUE-TRIAL EVENT RATE
# =============================================================================

actual_trial_event_rate = (
    len(actual_n_with_event) / TRIALS
)

control_trial_event_rate = (
    len(control_n_with_event) / TRIALS
)

print()
print("=" * 100)
print("TRIAL-LEVEL EVENT RATE")
print("=" * 100)

print(
    f"actual  = "
    f"{actual_trial_event_rate:.6f}"
)

print(
    f"control = "
    f"{control_trial_event_rate:.6f}"
)


# =============================================================================
# t DISTRIBUTION
# =============================================================================

print()
print("=" * 100)
print("t DISTRIBUTION: ACTUAL vs CONTROL")
print("=" * 100)

for t in range(T_MIN, T_MAX + 1):

    if t == 0:
        continue

    a = actual_t[t]
    c = control_t[t]

    if a == 0 and c == 0:
        continue

    print(
        f"t={t:+3} "
        f"actual={a:4} "
        f"control={c:4} "
        f"delta={a-c:+4}"
    )


# =============================================================================
# SYMMETRIC COMPARISON
# =============================================================================

print()
print("=" * 100)
print("SYMMETRIC |t| COMPARISON")
print("=" * 100)

for magnitude in range(
    1,
    max(abs(T_MIN), abs(T_MAX)) + 1
):

    an = actual_t[-magnitude]
    ap = actual_t[+magnitude]

    cn = control_t[-magnitude]
    cp = control_t[+magnitude]

    if (
        an == 0
        and ap == 0
        and cn == 0
        and cp == 0
    ):
        continue

    print(
        f"|t|={magnitude:3} "
        f"actual(-)={an:3} "
        f"actual(+)={ap:3} "
        f"control(-)={cn:3} "
        f"control(+)={cp:3}"
    )


# =============================================================================
# ACTUAL FACTOR-SHIFT GEOMETRY
# =============================================================================

print()
print("=" * 100)
print("ACTUAL FACTOR-SHIFT GEOMETRY")
print("=" * 100)

geometry = {}

for event in actual_events:

    dp = event["a"] - event["p"]
    dq = event["b"] - event["q"]

    geometry.setdefault(
        event["t"],
        []
    ).append(
        (dp, dq)
    )

for t in sorted(geometry):

    records = geometry[t]

    dp_values = [
        dp
        for dp, dq in records
    ]

    dq_values = [
        dq
        for dp, dq in records
    ]

    print()
    print(
        f"t={t:+3} "
        f"count={len(records)}"
    )

    print(
        f"    dp = "
        f"{min(dp_values):+,} ... "
        f"{max(dp_values):+,}"
    )

    print(
        f"    dq = "
        f"{min(dq_values):+,} ... "
        f"{max(dq_values):+,}"
    )


# =============================================================================
# PRODUCT IDENTITY
# =============================================================================

print()
print("=" * 100)
print("IDENTITY CHECK")
print("=" * 100)

failures = 0

for event in actual_events:

    lhs = (
        event["a"] * event["b"]
        - event["p"] * event["q"]
    )

    rhs = STEP * event["t"]

    if lhs != rhs:

        failures += 1

print(
    f"actual identity failures = {failures}"
)


# =============================================================================
# CONTROL RESIDUE CHECK
# =============================================================================

print()
print("=" * 100)
print("CONTROL RESIDUE CONDITION")
print("=" * 100)

bad_controls = 0

for event in control_events:

    if math.gcd(event["n"], STEP) != 1:

        bad_controls += 1

print(
    f"control gcd(n,2M) failures = "
    f"{bad_controls}"
)


# =============================================================================
# NORMALIZED COMPARISON
# =============================================================================

print()
print("=" * 100)
print("NORMALIZED t FREQUENCIES")
print("=" * 100)

total_actual = len(actual_events)
total_control = len(control_events)

for t in range(T_MIN, T_MAX + 1):

    if t == 0:
        continue

    a = actual_t[t]
    c = control_t[t]

    af = (
        a / total_actual
        if total_actual
        else 0.0
    )

    cf = (
        c / total_control
        if total_control
        else 0.0
    )

    if a or c:

        print(
            f"t={t:+3} "
            f"actual_freq={af:.8f} "
            f"control_freq={cf:.8f}"
        )


# =============================================================================
# EXTREME ACTUAL EVENTS
# =============================================================================

print()
print("=" * 100)
print("EXTREME ACTUAL EVENTS")
print("=" * 100)

if actual_events:

    by_abs_t = sorted(
        actual_events,
        key=lambda e: abs(e["t"]),
        reverse=True
    )[:20]

    print()
    print("Largest |t|:")

    for e in by_abs_t:

        dp = e["a"] - e["p"]
        dq = e["b"] - e["q"]

        print(
            f"t={e['t']:+3} "
            f"actual=({e['p']:,},{e['q']:,}) "
            f"shifted=({e['a']:,},{e['b']:,}) "
            f"dp={dp:+,} "
            f"dq={dq:+,}"
        )


# =============================================================================
# FINAL INTERPRETATION DATA
# =============================================================================

print()
print("=" * 100)
print("INTERPRETATION")
print("=" * 100)

print(
    "The actual experiment samples:"
)

print(
    "    N_t = p*q + 2*M*t"
)

print()

print(
    "The control samples the same t-lattice shape"
)

print(
    "    N_t = n0 + 2*M*t"
)

print()

print(
    "with gcd(n0,2M)=1."
)

print()

print(
    "Compare actual/control event rates before "
    "drawing conclusions about special collision structure."
)
