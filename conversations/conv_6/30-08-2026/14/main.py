#!/usr/bin/env python3
import math
import random
from itertools import combinations
from collections import Counter, defaultdict

# ============================================================================
# CONFIGURATION
# ============================================================================
M = 111_546_435
FACTOR_MIN = 10_000
FACTOR_MAX = 100_000
TRIALS = 300
SEED = 1_511_464_998

# Candidate modulus primes.  We search triples r1<r2<r3 with R=r1*r2*r3<n.
MOD_MIN = 300
MOD_MAX = 3_000

# Keep the search practical.
CLOSE_RATIO = 0.20       # each modulus must be within 20% of the smallest
RANDOM_TRIPLES = 20      # per anchor for matched-random control
TOP_EXAMPLES = 15

# ============================================================================
# PRIME UTILITIES
# ============================================================================
def sieve(limit: int):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, int(limit ** 0.5) + 1):
        if sieve[p]:
            sieve[p*p:limit+1:p] = b"\x00" * (((limit - p*p)//p) + 1)
    return [i for i, v in enumerate(sieve) if v]


def is_prime_small(n: int, primes):
    return n in primes

# ============================================================================
# ANCHORS
# ============================================================================
def make_anchors(prime_pool, count, rng):
    anchors = []
    seen = set()
    while len(anchors) < count:
        p, q = rng.sample(prime_pool, 2)
        if p > q:
            p, q = q, p
        key = (p, q)
        if key in seen:
            continue
        seen.add(key)
        anchors.append((p, q, p*q))
    return anchors

# ============================================================================
# MODULUS TRIPLES
# ============================================================================
def build_modulus_primes():
    all_primes = sieve(MOD_MAX)
    return [p for p in all_primes if MOD_MIN <= p <= MOD_MAX]


def best_max_product_triple(n, mod_primes):
    # Find max R<n among close triples.  Since R must be close to n,
    # start around cube-root(n) and only scan the nearby tail.
    target = n ** (1.0 / 3.0)
    idx = 0
    while idx < len(mod_primes) and mod_primes[idx] < target * (1.0 - CLOSE_RATIO):
        idx += 1

    best = None
    # A modest local window is enough because the objective is product-first.
    start = max(0, idx - 25)
    end = min(len(mod_primes), idx + 90)
    pool = mod_primes[start:end]

    # For correctness, fall back to a somewhat wider scan if needed.
    if len(pool) < 10:
        pool = mod_primes

    for i in range(len(pool)):
        r1 = pool[i]
        for j in range(i + 1, len(pool)):
            r2 = pool[j]
            if r2 > int(r1 * (1.0 + CLOSE_RATIO)):
                break
            r12 = r1 * r2
            if r12 >= n:
                break
            # Need r3 >= r2 and close to r1.
            for k in range(j + 1, len(pool)):
                r3 = pool[k]
                if r3 > int(r1 * (1.0 + CLOSE_RATIO)):
                    break
                R = r12 * r3
                if R >= n:
                    break
                if best is None or R > best[0]:
                    best = (R, (r1, r2, r3))
    return best


def random_valid_triple(rng, mod_primes, n):
    # Generate random close triple satisfying R<n.
    for _ in range(200):
        r1 = rng.choice(mod_primes)
        lo = r1
        hi = min(MOD_MAX, int(r1 * (1.0 + CLOSE_RATIO)))
        candidates = [r for r in mod_primes if lo < r <= hi]
        if len(candidates) < 2:
            continue
        r2, r3 = sorted(rng.sample(candidates, 2))
        R = r1 * r2 * r3
        if R < n:
            return R, tuple(sorted((r1, r2, r3)))
    return None

# ============================================================================
# FINGERPRINT / BRANCH LOGIC
# ============================================================================
def branch_from_factor_pair(p, q):
    a, b = p % 4, q % 4
    if a > b:
        a, b = b, a
    if (a, b) == (1, 1):
        return "11"
    if (a, b) == (1, 3):
        return "13/31"
    if (a, b) == (3, 3):
        return "33"
    raise ValueError((p, q, a, b))


def allowed_branches_for_n(n):
    # Odd n can only be 1 or 3 mod 4.
    if n % 4 == 1:
        return ("11", "33")
    if n % 4 == 3:
        return ("13/31",)
    raise ValueError("n must be odd")


def pair_product_ok(a, b, n, mods):
    return all((a * b) % r == n % r for r in mods)


def pair_candidates(n, mods, prime_pool, branch_filter=None):
    # Enumerate p only.  For each p, each modulus constrains q's residue.
    # We simply test every q prime using modular product checks, but exploit
    # a cheap first-modulus residue filter before doing all moduli.
    r0 = mods[0]
    target0 = n % r0
    by_res = defaultdict(list)
    for q in prime_pool:
        by_res[q % r0].append(q)

    out = []
    for a in prime_pool:
        if branch_filter is not None and branch_from_factor_pair(a, a) not in branch_filter:
            # Only a cheap parity-class precheck.  We still need q's branch.
            pass
        ra = a % r0
        candidates_q = by_res.get((target0 * pow(ra, -1, r0)) % r0, ()) if math.gcd(ra, r0) == 1 else ()
        for b in candidates_q:
            if a > b:
                continue
            if branch_filter is not None and branch_from_factor_pair(a, b) not in branch_filter:
                continue
            if pair_product_ok(a, b, n, mods):
                out.append((a, b))
    return out


def p_fingerprint_collisions(p, prime_pool, mods):
    fp = tuple(p % r for r in mods)
    matches = [q for q in prime_pool if tuple(q % r for r in mods) == fp]
    return matches

# ============================================================================
# MAIN EXPERIMENT
# ============================================================================
def run():
    print("=" * 100)
    print("THREE-MODULUS C-FIRST / P-FINGERPRINT / FACTOR-PAIR RESOLUTION EXPERIMENT")
    print("=" * 100)
    print(f"M                         = {M:,}")
    print(f"factor range             = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"actual anchors            = {TRIALS}")
    print(f"modulus prime range      = {MOD_MIN} - {MOD_MAX}")
    print(f"close ratio              = {CLOSE_RATIO:.0%}")
    print(f"random triples / anchor  = {RANDOM_TRIPLES}")
    print(f"seed                      = {SEED:,}")
    print()

    rng = random.Random(SEED)
    primes_all = sieve(FACTOR_MAX)
    prime_pool = [p for p in primes_all if FACTOR_MIN <= p <= FACTOR_MAX]
    mod_primes = build_modulus_primes()

    print("=" * 100)
    print("PRIME POOLS")
    print("=" * 100)
    print(f"factor primes             = {len(prime_pool):,}")
    print(f"modulus primes            = {len(mod_primes):,}")

    anchors = make_anchors(prime_pool, TRIALS, rng)
    print(f"actual anchors             = {len(anchors)}")

    # ------------------------------------------------------------------------
    # C: choose r1,r2,r3 and test mod-4 branch first.
    # ------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("C-FIRST: MOD-4 BRANCH DISCRIMINATION")
    print("=" * 100)

    records = []
    branch_counts = Counter()
    branch_resolved = 0

    for idx, (p, q, n) in enumerate(anchors, 1):
        best = best_max_product_triple(n, mod_primes)
        if best is None:
            raise RuntimeError(f"No valid triple for n={n}")
        R, mods = best

        allowed = allowed_branches_for_n(n)
        actual_branch = branch_from_factor_pair(p, q)

        # First, branch-only: enumerate compatible pairs in each branch.
        candidates_all = pair_candidates(n, mods, prime_pool, branch_filter=None)
        surviving = sorted(set(branch_from_factor_pair(a, b) for a, b in candidates_all))
        actual_survives = actual_branch in surviving
        if len(surviving) == 1:
            branch_resolved += 1

        rec = {
            "p": p, "q": q, "n": n,
            "R": R, "mods": mods,
            "R_over_n": R / n,
            "gap": n - R,
            "allowed_branches": allowed,
            "actual_branch": actual_branch,
            "surviving_branches": tuple(surviving),
            "all_pairs": tuple(candidates_all),
        }
        records.append(rec)
        branch_counts[actual_branch] += 1

        if idx % 50 == 0:
            print(f"anchor {idx:3d}/{len(anchors)}")

    print("\nBRANCH SUMMARY")
    print("-" * 100)
    print(f"anchors with exactly one surviving branch = {branch_resolved}/{TRIALS}")
    for b in ("11", "13/31", "33"):
        print(f"actual branch {b:5s} = {branch_counts[b]:3d}")

    print("\nBRANCH SURVIVAL COUNTS")
    print("-" * 100)
    print(" surviving branches     anchors")
    surv_counter = Counter(tuple(r["surviving_branches"]) for r in records)
    for k, v in sorted(surv_counter.items(), key=lambda kv: (len(kv[0]), kv[0])):
        print(f" {str(k):22s} {v:7d}")

    # ------------------------------------------------------------------------
    # A: unique p fingerprint.
    # ------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("A: UNIQUE p-FINGERPRINT")
    print("=" * 100)
    print("For the selected (r1,r2,r3): fingerprint(p)=(p mod r1,p mod r2,p mod r3).")
    print()

    p_unique = 0
    p_ambig = 0
    p_max_bucket = 0
    p_bucket_sizes = Counter()

    for rec in records:
        matches = p_fingerprint_collisions(rec["p"], prime_pool, rec["mods"])
        rec["p_matches"] = tuple(matches)
        k = len(matches)
        p_bucket_sizes[k] += 1
        p_max_bucket = max(p_max_bucket, k)
        if k == 1:
            p_unique += 1
        else:
            p_ambig += 1

    print(f"unique p fingerprints      = {p_unique}/{TRIALS}")
    print(f"ambiguous p fingerprints   = {p_ambig}/{TRIALS}")
    print(f"maximum collision bucket   = {p_max_bucket}")
    print("bucket distribution:")
    for k in sorted(p_bucket_sizes):
        print(f"  {k:3d} primes sharing fingerprint : {p_bucket_sizes[k]:3d} anchors")

    # ------------------------------------------------------------------------
    # B: factor-pair recovery from n residues.
    # ------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("B: FACTOR-PAIR RECOVERY FROM n MOD (r1,r2,r3)")
    print("=" * 100)

    pair_unique = 0
    pair_ambig = 0
    pair_zero = 0
    cand_hist = Counter()

    for rec in records:
        cands = rec["all_pairs"]
        cand_hist[len(cands)] += 1
        if len(cands) == 1 and cands[0] == (rec["p"], rec["q"]):
            pair_unique += 1
        elif len(cands) == 0:
            pair_zero += 1
        else:
            pair_ambig += 1

    print(f"unique actual factor pair   = {pair_unique}/{TRIALS}")
    print(f"ambiguous factor pair       = {pair_ambig}/{TRIALS}")
    print(f"zero compatible pairs       = {pair_zero}/{TRIALS}")
    print("candidate-count distribution:")
    for k in sorted(cand_hist):
        print(f"  {k:3d} candidates : {cand_hist[k]:3d} anchors")

    # ------------------------------------------------------------------------
    # Combined diagnostic: C -> A -> B
    # ------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("C -> A -> B COMBINED RESULT")
    print("=" * 100)
    print("C = branch discrimination")
    print("A = unique p fingerprint")
    print("B = unique factor pair")
    print()
    print(f"C unique branch        = {branch_resolved:3d}/{TRIALS}")
    print(f"A unique p             = {p_unique:3d}/{TRIALS}")
    print(f"B unique factor pair   = {pair_unique:3d}/{TRIALS}")

    # ------------------------------------------------------------------------
    # Continue: compare max-product triples against matched random triples.
    # ------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("CONTINUATION: MAX-PRODUCT vs RANDOM CLOSE TRIPLES")
    print("=" * 100)

    actual_ratios = []
    random_ratios = []
    actual_unique = 0
    actual_ambig = 0
    random_unique = 0
    random_ambig = 0
    matched_rows = []

    for i, rec in enumerate(records, 1):
        if len(rec["all_pairs"]) == 1:
            actual_unique += 1
        else:
            actual_ambig += 1
        actual_ratios.append(rec["R_over_n"])

        # Only sample random triples whose R/n is reasonably comparable.
        accepted = 0
        for _ in range(RANDOM_TRIPLES * 10):
            got = random_valid_triple(rng, mod_primes, rec["n"])
            if got is None:
                continue
            Rr, mr = got
            ratio = Rr / rec["n"]
            if ratio >= max(0.0, rec["R_over_n"] - 0.10):
                cands = pair_candidates(rec["n"], mr, prime_pool, branch_filter=None)
                random_ratios.append(ratio)
                if len(cands) == 1:
                    random_unique += 1
                else:
                    random_ambig += 1
                accepted += 1
                if accepted >= RANDOM_TRIPLES:
                    break
        if i % 50 == 0:
            print(f"anchor {i:3d}/{len(records)}")

    def mean(xs):
        return sum(xs) / len(xs) if xs else float("nan")

    print("\nSUMMARY")
    print("-" * 100)
    print(f"actual triples                 = {len(actual_ratios):,}")
    print(f"random accepted triples        = {len(random_ratios):,}")
    print(f"actual mean R/n                = {mean(actual_ratios):.9f}")
    print(f"random mean R/n                = {mean(random_ratios):.9f}")
    print(f"actual unique factor pairs     = {actual_unique}/{len(actual_ratios)} ({actual_unique/max(1,len(actual_ratios)):.2%})")
    print(f"random unique factor pairs     = {random_unique}/{len(random_ratios)} ({random_unique/max(1,len(random_ratios)):.2%})")

    # ------------------------------------------------------------------------
    # Top examples.
    # ------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("TOP MAX-PRODUCT EXAMPLES")
    print("=" * 100)
    for rec in sorted(records, key=lambda r: r["R_over_n"], reverse=True)[:TOP_EXAMPLES]:
        print(
            f"n={rec['n']:,} p={rec['p']:,} q={rec['q']:,} "
            f"mods={rec['mods']} R={rec['R']:,} R/n={rec['R_over_n']:.10f} "
            f"gap={rec['gap']:,} branch={rec['actual_branch']} "
            f"branches={rec['surviving_branches']} pBucket={len(rec['p_matches'])} "
            f"pairCandidates={len(rec['all_pairs'])}"
        )

    print("\n" + "=" * 100)
    print("IMPORTANT INTERPRETATION")
    print("=" * 100)
    print("1. C is tested first: does n mod 4 plus the three product residues eliminate a branch?")
    print("2. A asks whether the three residues uniquely identify p among primes in the factor interval.")
    print("3. B asks the stronger question: can n mod r1,r2,r3 uniquely identify the whole factor pair?")
    print("4. The random control checks whether any advantage is simply caused by R being close to n.")
    print()
    print("Note: uniqueness is conditional on the tested prime interval and the selected modulus triple.")
    print("It is not a general statement that three residues factor arbitrary integers.")
    print()
    print("Experiment complete.")

if __name__ == "__main__":
    run()
