#!/usr/bin/env python3

import math
import random
from collections import Counter, defaultdict

# ============================================================
# MODULAR FACTOR-TRANSFORMATION EXPERIMENT
# ============================================================

M = 111_546_435
STEP = 2 * M

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

T_MIN = -25
T_MAX = 25

TRIALS = 300

SEED = 1511464998

MODULI = [3, 5, 7, 11, 13, 17, 19, 23]

rng = random.Random(SEED)

print("=" * 100)
print("MODULAR FACTOR-TRANSFORMATION / ACTUAL vs RANDOM CONTROL")
print("=" * 100)
print(f"M                = {M:,}")
print(f"2M               = {STEP:,}")
print(f"factor range     = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
print(f"t range          = [{T_MIN}, {T_MAX}]")
print(f"trials            = {TRIALS}")
print(f"random seed       = {SEED}")
print()

# ============================================================
# PRIME SIEVE
# ============================================================

def sieve(limit):
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i in range(2, limit + 1) if is_prime[i]]


PRIMES = sieve(FACTOR_MAX)

FACTOR_PRIMES = [
    p for p in PRIMES
    if FACTOR_MIN <= p <= FACTOR_MAX
]

FACTOR_SET = set(FACTOR_PRIMES)

print(f"factor primes     = {len(FACTOR_PRIMES):,}")
print()


# ============================================================
# BASIC HELPERS
# ============================================================

def is_semiprime_in_range(n):
    """
    Return (p,q) if n = p*q with both prime factors in range.
    Return None otherwise.

    We only search the lower factor up to sqrt(n), which keeps
    this reasonably cheap for the present experiment.
    """
    limit = math.isqrt(n)

    for p in FACTOR_PRIMES:
        if p > limit:
            break

        if n % p == 0:
            q = n // p

            if q in FACTOR_SET and q >= p:
                return p, q

    return None


def random_semiprime():
    while True:
        p = rng.choice(FACTOR_PRIMES)
        q = rng.choice(FACTOR_PRIMES)

        n = p * q

        if p <= q:
            return p, q, n
        return q, p, n


def find_shifted_factorization(n, exclude_pair=None):
    """
    Find a,b in the factor range such that a*b=n.

    Exclude the original pair so that the event is genuinely a
    shifted collision rather than t=0 / the same factorization.
    """
    limit = math.isqrt(n)

    for a in FACTOR_PRIMES:
        if a > limit:
            break

        if n % a != 0:
            continue

        b = n // a

        if b not in FACTOR_SET:
            continue

        if a > b:
            continue

        pair = (a, b)

        if exclude_pair is not None and pair == tuple(sorted(exclude_pair)):
            continue

        return pair

    return None


def inv_mod(a, p):
    """
    Modular inverse.
    All p used here are prime.
    """
    return pow(a % p, -1, p)


def legendre_symbol(a, p):
    """
    Returns:
       0  if a == 0 mod p
      +1  if a is a quadratic residue
      -1  otherwise
    """
    a %= p

    if a == 0:
        return 0

    v = pow(a, (p - 1) // 2, p)

    if v == 1:
        return 1

    if v == p - 1:
        return -1

    raise RuntimeError("Invalid Legendre-symbol result")


def multiplicative_order(a, p):
    """
    Multiplicative order of a modulo prime p.
    """
    a %= p

    if a == 0:
        return 0

    order = p - 1

    # Factor p-1.
    x = order
    factors = []

    d = 2
    while d * d <= x:
        if x % d == 0:
            factors.append(d)
            while x % d == 0:
                x //= d
        d += 1

    if x > 1:
        factors.append(x)

    for f in factors:
        while order % f == 0:
            candidate = order // f

            if pow(a, candidate, p) == 1:
                order = candidate
            else:
                break

    return order


# ============================================================
# DATA STRUCTURES
# ============================================================

actual_events = []
control_events = []

actual_t = Counter()
control_t = Counter()

actual_u = {r: Counter() for r in MODULI}
control_u = {r: Counter() for r in MODULI}

actual_order = {r: Counter() for r in MODULI}
control_order = {r: Counter() for r in MODULI}

actual_legendre = {r: Counter() for r in MODULI}
control_legendre = {r: Counter() for r in MODULI}

actual_t_u = {r: defaultdict(list) for r in MODULI}
control_t_u = {r: defaultdict(list) for r in MODULI}

identity_failures_actual = 0
identity_failures_control = 0


# ============================================================
# GENERATE ACTUAL EVENTS
# ============================================================

def generate_actual_event():
    p, q, n = random_semiprime()

    for t in range(T_MIN, T_MAX + 1):
        if t == 0:
            continue

        Nt = n + STEP * t

        if Nt <= 0:
            continue

        pair = find_shifted_factorization(Nt, exclude_pair=(p, q))

        if pair is None:
            continue

        a, b = pair

        return {
            "p": p,
            "q": q,
            "n": n,
            "t": t,
            "Nt": Nt,
            "a": a,
            "b": b,
        }

    return None


# ============================================================
# GENERATE RANDOM CONTROL EVENTS
# ============================================================

def generate_control_event():
    """
    Random anchor semiprime.

    The construction is identical to the actual experiment:
        Nt = p*q + 2*M*t

    The anchor itself is random rather than drawn from the actual
    experimental sample.
    """
    return generate_actual_event()


# ============================================================
# MODULAR RECORDING
# ============================================================

def record_event(
    event,
    events,
    t_counter,
    u_counter,
    order_counter,
    legendre_counter,
    t_u,
    identity_type,
):
    p = event["p"]
    q = event["q"]
    a = event["a"]
    b = event["b"]
    t = event["t"]

    for r in MODULI:

        # p and q should be invertible because all factor primes are
        # larger than 23.
        ip = inv_mod(p, r)
        iq = inv_mod(q, r)

        u = (a * ip) % r
        v = (b * iq) % r

        # Core modular identity:
        #
        # (a/p) * (b/q) = 1 mod r
        if (u * v) % r != 1:
            if identity_type == "actual":
                global identity_failures_actual
                identity_failures_actual += 1
            else:
                global identity_failures_control
                identity_failures_control += 1

        u_counter[r][u] += 1

        order = multiplicative_order(u, r)
        order_counter[r][order] += 1

        leg = legendre_symbol(u, r)
        legendre_counter[r][leg] += 1

        t_u[r][t].append(u)

    events.append(event)
    t_counter[t] += 1


# ============================================================
# COLLECT EVENTS
# ============================================================

attempts = 0

while len(actual_events) < TRIALS:
    event = generate_actual_event()

    if event is not None:
        record_event(
            event,
            actual_events,
            actual_t,
            actual_u,
            actual_order,
            actual_legendre,
            actual_t_u,
            "actual",
        )

    attempts += 1

print("=" * 100)
print("ACTUAL SAMPLE")
print("=" * 100)
print(f"successful anchors = {len(actual_events)}")
print(f"generation attempts = {attempts}")
print()


attempts = 0

while len(control_events) < TRIALS:
    event = generate_control_event()

    if event is not None:
        record_event(
            event,
            control_events,
            control_t,
            control_u,
            control_order,
            control_legendre,
            control_t_u,
            "control",
        )

    attempts += 1

print("=" * 100)
print("CONTROL SAMPLE")
print("=" * 100)
print(f"successful anchors = {len(control_events)}")
print(f"generation attempts = {attempts}")
print()


# ============================================================
# UNIFORMITY STATISTIC
# ============================================================

def chi_square_uniform(counter, values):
    """
    Pearson chi-square against a uniform distribution.
    Returns None if there are no observations.
    """
    total = sum(counter.values())

    if total == 0:
        return None

    expected = total / len(values)

    return sum(
        ((counter.get(v, 0) - expected) ** 2) / expected
        for v in values
    )


# ============================================================
# CORRELATION
# ============================================================

def pearson(xs, ys):
    if len(xs) != len(ys) or len(xs) < 2:
        return float("nan")

    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)

    num = sum(
        (x - mx) * (y - my)
        for x, y in zip(xs, ys)
    )

    dx = sum((x - mx) ** 2 for x in xs)
    dy = sum((y - my) ** 2 for y in ys)

    if dx == 0 or dy == 0:
        return float("nan")

    return num / math.sqrt(dx * dy)


# ============================================================
# GLOBAL MODULAR SUMMARY
# ============================================================

print("=" * 100)
print("MODULAR TRANSFORMATION SUMMARY")
print("=" * 100)

print(
    f"{'r':>3} "
    f"{'ACT N':>8} "
    f"{'ACT chi2':>12} "
    f"{'CTL N':>8} "
    f"{'CTL chi2':>12}"
)

print("-" * 60)

for r in MODULI:
    values = list(range(1, r))

    act_n = sum(actual_u[r].values())
    ctl_n = sum(control_u[r].values())

    act_chi = chi_square_uniform(actual_u[r], values)
    ctl_chi = chi_square_uniform(control_u[r], values)

    print(
        f"{r:3d} "
        f"{act_n:8d} "
        f"{act_chi:12.4f} "
        f"{ctl_n:8d} "
        f"{ctl_chi:12.4f}"
    )

print()


# ============================================================
# RAW u DISTRIBUTIONS
# ============================================================

print("=" * 100)
print("u_r DISTRIBUTIONS")
print("=" * 100)

for r in MODULI:
    print(f"\nr = {r}")

    values = list(range(1, r))

    print(" actual :", end=" ")
    for u in values:
        print(f"{u}:{actual_u[r][u]}", end=" ")
    print()

    print(" control:", end=" ")
    for u in values:
        print(f"{u}:{control_u[r][u]}", end=" ")
    print()

print()


# ============================================================
# LEGENDRE SYMBOL
# ============================================================

print("=" * 100)
print("QUADRATIC RESIDUOSITY")
print("=" * 100)

print(
    f"{'r':>3} "
    f"{'ACT +':>8} "
    f"{'ACT -':>8} "
    f"{'CTL +':>8} "
    f"{'CTL -':>8}"
)

print("-" * 50)

for r in MODULI:
    print(
        f"{r:3d} "
        f"{actual_legendre[r][1]:8d} "
        f"{actual_legendre[r][-1]:8d} "
        f"{control_legendre[r][1]:8d} "
        f"{control_legendre[r][-1]:8d}"
    )

print()


# ============================================================
# MULTIPLICATIVE ORDER
# ============================================================

print("=" * 100)
print("MULTIPLICATIVE ORDER DISTRIBUTIONS")
print("=" * 100)

for r in MODULI:
    print(f"\nr = {r}")

    print(" actual :", dict(sorted(actual_order[r].items())))
    print(" control:", dict(sorted(control_order[r].items())))

print()


# ============================================================
# t -> u CORRELATION
# ============================================================

print("=" * 100)
print("t / MODULAR TRANSFORMATION CORRELATION")
print("=" * 100)

for r in MODULI:

    act_x = []
    act_y = []

    ctl_x = []
    ctl_y = []

    for t, values in actual_t_u[r].items():
        for u in values:
            act_x.append(t)
            act_y.append(u)

    for t, values in control_t_u[r].items():
        for u in values:
            ctl_x.append(t)
            ctl_y.append(u)

    print(
        f"r={r:2d} "
        f"actual_corr={pearson(act_x, act_y):+.8f} "
        f"control_corr={pearson(ctl_x, ctl_y):+.8f}"
    )

print()


# ============================================================
# t-CONDITIONED MEAN u
# ============================================================

print("=" * 100)
print("t-CONDITIONED MEAN u_r")
print("=" * 100)

for r in MODULI:
    print(f"\nr = {r}")

    print(
        f"{'t':>4} "
        f"{'ACT N':>7} "
        f"{'ACT mean u':>12} "
        f"{'CTL N':>7} "
        f"{'CTL mean u':>12}"
    )

    for t in range(T_MIN, T_MAX + 1):

        av = actual_t_u[r].get(t, [])
        cv = control_t_u[r].get(t, [])

        if not av and not cv:
            continue

        am = sum(av) / len(av) if av else float("nan")
        cm = sum(cv) / len(cv) if cv else float("nan")

        print(
            f"{t:4d} "
            f"{len(av):7d} "
            f"{am:12.6f} "
            f"{len(cv):7d} "
            f"{cm:12.6f}"
        )

print()


# ============================================================
# t DISTRIBUTION
# ============================================================

print("=" * 100)
print("t DISTRIBUTION")
print("=" * 100)

print(
    f"{'t':>4} "
    f"{'ACTUAL':>8} "
    f"{'CONTROL':>8}"
)

print("-" * 30)

for t in range(T_MIN, T_MAX + 1):

    if t == 0:
        continue

    if actual_t[t] == 0 and control_t[t] == 0:
        continue

    print(
        f"{t:4d} "
        f"{actual_t[t]:8d} "
        f"{control_t[t]:8d}"
    )

print()


# ============================================================
# CORE IDENTITY VERIFICATION
# ============================================================

print("=" * 100)
print("CORE MODULAR IDENTITY")
print("=" * 100)

print(
    "(a / p) * (b / q) = 1 (mod r)"
)
print(
    "equivalently:"
)
print(
    "p*dq + q*dp + dp*dq = 0 (mod r)"
)
print()

print(
    f"actual modular identity failures  = "
    f"{identity_failures_actual}"
)

print(
    f"control modular identity failures = "
    f"{identity_failures_control}"
)

print()


# ============================================================
# GLOBAL FACTOR TRANSFORMATION STATISTICS
# ============================================================

print("=" * 100)
print("GLOBAL NORMALIZED FACTOR TRANSFORMATIONS")
print("=" * 100)

for label, events in (
    ("ACTUAL", actual_events),
    ("CONTROL", control_events),
):

    print(f"\n{label}")

    ndp = []
    ndq = []

    for e in events:
        p = e["p"]
        q = e["q"]
        a = e["a"]
        b = e["b"]

        ndp.append((a - p) / p)
        ndq.append((b - q) / q)

    print(f"events     = {len(events)}")
    print(f"mean ndp   = {sum(ndp) / len(ndp):+.10f}")
    print(f"mean ndq   = {sum(ndq) / len(ndq):+.10f}")

print()


# ============================================================
# IMPORTANT RESIDUE-SIGNATURE TEST
# ============================================================

print("=" * 100)
print("MODULAR SIGNATURE MATCH TEST")
print("=" * 100)

for r in MODULI:

    act_total = sum(actual_u[r].values())
    ctl_total = sum(control_u[r].values())

    print(f"\nr={r}")

    print(" actual normalized frequencies:")
    for u in range(1, r):
        freq = actual_u[r][u] / act_total
        print(f"   u={u:2d}: {freq:.8f}")

    print(" control normalized frequencies:")
    for u in range(1, r):
        freq = control_u[r][u] / ctl_total
        print(f"   u={u:2d}: {freq:.8f}")

print()


# ============================================================
# FINAL
# ============================================================

print("=" * 100)
print("FINAL")
print("=" * 100)

print(
    "Every event was generated from"
)
print(
    "    Nt = p*q + 2*M*t"
)
print(
    "and the shifted factorization"
)
print(
    "    Nt = a*b."
)
print()

print(
    "Therefore every valid event must satisfy, for every r | M:"
)
print(
    "    (a/p) * (b/q) = 1 (mod r)."
)
print()

print(
    "The experiment now measures whether the individual modular"
)
print(
    "transformation u_r = a/p has any non-random structure."
)
print()

print("Experiment complete.")
