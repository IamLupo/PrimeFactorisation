#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include <gmpxx.h>

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct PrimePair {
    u64 p;
    u64 q;
    u64 n;
    u64 s;
};

struct Witness {
    bool exists = false;
    u64 m = 0;
    u64 k = 0;
    u64 t = 0;
    int sign = 0; // -1 => kt = mr-1, +1 => kt = mr+1
};

struct Event {
    u64 K = 0;
    u64 m = 0;
    u64 k = 0;
    u64 t = 0;
    int sign = 0;
};

struct Stats {
    u64 cases = 0;

    u64 winner_changes = 0;

    u64 changes_on_new_divisor = 0;
    u64 changes_without_new_divisor = 0;

    u64 winner_stays_between_events = 0;

    u64 record_events = 0;

    u64 m1_records = 0;
    u64 m2_records = 0;
    u64 m3_records = 0;
    u64 m4_records = 0;

    u64 transition_1_to_2 = 0;
    u64 transition_2_to_3 = 0;
    u64 transition_3_to_4 = 0;

    u64 profile_failures = 0;
    u64 verification_failures = 0;

    u64 min_event_gap =
        std::numeric_limits<u64>::max();

    u64 max_event_gap = 0;

    u64 sum_event_gap = 0;

    u64 examples_printed = 0;
};

static mpz_class mpz_from_u64(u64 x) {
    return mpz_class(std::to_string(x));
}

static bool is_prime(u64 n) {
    if (n < 2) {
        return false;
    }

    const mpz_class z = mpz_from_u64(n);

    return mpz_probab_prime_p(
        z.get_mpz_t(),
        25
    ) > 0;
}

static u64 integer_sqrt(u64 n) {
    u64 lo = 0;
    u64 hi = 1000000000ULL;

    while (lo <= hi) {
        const u64 mid = lo + (hi - lo) / 2;

        const u128 sq =
            static_cast<u128>(mid) * mid;

        if (sq <= n) {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }

    return hi;
}

static u64 random_prime(
    std::mt19937_64& rng,
    u64 lo,
    u64 hi
) {
    std::uniform_int_distribution<u64> dist(lo, hi);

    while (true) {
        const u64 x = dist(rng) | 1ULL;

        if (x < lo || x > hi) {
            continue;
        }

        if (is_prime(x)) {
            return x;
        }
    }
}

static PrimePair generate_case(
    std::mt19937_64& rng
) {
    const u64 a =
        random_prime(rng, 10000, 1000000);

    u64 b =
        random_prime(rng, 10000, 1000000);

    while (b == a) {
        b =
            random_prime(rng, 10000, 1000000);
    }

    const u64 p = std::min(a, b);
    const u64 q = std::max(a, b);

    const u64 n = p * q;
    const u64 s = integer_sqrt(n);

    return {p, q, n, s};
}

/*
 * Compute the best witness for fixed m and K.

 * The largest divisor k <= K gives the smallest t.
 */
static Witness best_for_m(
    u64 r,
    u64 K,
    u64 m
) {
    Witness best;

    const u128 mr =
        static_cast<u128>(m) * r;

    for (int sign_index = 0;
         sign_index < 2;
         ++sign_index) {

        const int sign =
            sign_index == 0 ? -1 : +1;

        u128 value = 0;

        if (sign == -1) {
            if (mr <= 1) {
                continue;
            }

            value = mr - 1;
        } else {
            value = mr + 1;
        }

        if (
            value >
            static_cast<u128>(
                std::numeric_limits<u64>::max()
            )
        ) {
            continue;
        }

        const u64 value_u64 =
            static_cast<u64>(value);

        for (u64 k = K; k >= 2; --k) {
            if (value_u64 % k != 0) {
                continue;
            }

            const u64 t =
                value_u64 / k;

            if (t == 0) {
                continue;
            }

            if (t > r / 2) {
                continue;
            }

            if (
                !best.exists ||
                t < best.t ||
                (t == best.t && k < best.k)
            ) {
                best.exists = true;
                best.m = m;
                best.k = k;
                best.t = t;
                best.sign = sign;
            }

            break;
        }
    }

    return best;
}

static Witness profile_winner(
    u64 r,
    u64 K
) {
    Witness best;

    for (u64 m = 1; m <= 4; ++m) {
        const Witness w =
            best_for_m(r, K, m);

        if (!w.exists) {
            continue;
        }

        if (
            !best.exists ||
            w.t < best.t ||
            (
                w.t == best.t &&
                w.m < best.m
            )
        ) {
            best = w;
        }
    }

    return best;
}

static bool verify_witness(
    u64 r,
    const Witness& w
) {
    if (!w.exists) {
        return false;
    }

    const u128 lhs =
        static_cast<u128>(w.k) * w.t;

    const u128 rhs =
        static_cast<u128>(w.m) * r;

    if (w.sign == +1) {
        return lhs == rhs + 1;
    }

    if (w.sign == -1) {
        return lhs + 1 == rhs;
    }

    return false;
}

/*
 * Check whether K introduces a new divisor for any m.

 * A new divisor event occurs exactly when K divides
 * one of:
 *
 *     m*r-1
 *     m*r+1
 *
 * for some m <= 4.
 */
static bool is_new_divisor_event(
    u64 r,
    u64 K
) {
    for (u64 m = 1; m <= 4; ++m) {
        const u128 mr =
            static_cast<u128>(m) * r;

        if (mr > K) {
            const u128 minus_value = mr - 1;

            if (
                minus_value %
                K == 0
            ) {
                return true;
            }
        }

        const u128 plus_value = mr + 1;

        if (
            plus_value %
            K == 0
        ) {
            return true;
        }
    }

    return false;
}

static void count_record(
    u64 old_k,
    const Witness& old_w,
    const Witness& new_w,
    Stats& stats
) {
    if (!new_w.exists) {
        return;
    }

    if (
        !old_w.exists ||
        new_w.k != old_w.k ||
        new_w.m != old_w.m ||
        new_w.t != old_w.t ||
        new_w.sign != old_w.sign
    ) {
        ++stats.record_events;

        switch (new_w.m) {
            case 1:
                ++stats.m1_records;
                break;

            case 2:
                ++stats.m2_records;
                break;

            case 3:
                ++stats.m3_records;
                break;

            case 4:
                ++stats.m4_records;
                break;

            default:
                break;
        }

        if (old_k != 0) {
            const u64 gap =
                new_w.k >= old_k
                    ? new_w.k - old_k
                    : 0;

            stats.min_event_gap =
                std::min(
                    stats.min_event_gap,
                    gap
                );

            stats.max_event_gap =
                std::max(
                    stats.max_event_gap,
                    gap
                );

            stats.sum_event_gap += gap;
        }
    }
}

static void analyze_prime(
    u64 r,
    u64 K_LIMIT,
    Stats& stats,
    bool print_examples
) {
    ++stats.cases;

    Witness previous;
    u64 previous_K = 0;

    bool previous_was_event = false;

    for (u64 K = 2;
         K <= K_LIMIT;
         ++K) {

        const Witness current =
            profile_winner(r, K);

        if (!current.exists) {
            continue;
        }

        if (!verify_witness(r, current)) {
            ++stats.verification_failures;
        }

        const bool event =
            is_new_divisor_event(
                r,
                K
            );

        if (event) {
            previous_was_event = true;
        }

        if (previous.exists) {
            const bool changed =
                (
                    current.m != previous.m ||
                    current.k != previous.k ||
                    current.t != previous.t ||
                    current.sign != previous.sign
                );

            if (changed) {
                ++stats.winner_changes;

                if (event) {
                    ++stats.changes_on_new_divisor;

                    if (
                        print_examples &&
                        stats.examples_printed < 15
                    ) {
                        ++stats.examples_printed;

                        std::cout
                            << "WINNER_CHANGE_ON_EVENT"
                            << " r=" << r
                            << " K=" << K
                            << " OLD_M=" << previous.m
                            << " OLD_K=" << previous.k
                            << " OLD_T=" << previous.t
                            << " NEW_M=" << current.m
                            << " NEW_K=" << current.k
                            << " NEW_T=" << current.t
                            << '\n';
                    }
                } else {
                    ++stats.changes_without_new_divisor;
                }

                if (
                    previous.m == 1 &&
                    current.m == 2
                ) {
                    ++stats.transition_1_to_2;
                }

                if (
                    previous.m == 2 &&
                    current.m == 3
                ) {
                    ++stats.transition_2_to_3;
                }

                if (
                    previous.m == 3 &&
                    current.m == 4
                ) {
                    ++stats.transition_3_to_4;
                }
            } else if (
                event &&
                previous_was_event
            ) {
                ++stats.winner_stays_between_events;
            }
        }

        count_record(
            previous_K,
            previous,
            current,
            stats
        );

        previous = current;
        previous_K = K;
        previous_was_event = event;
    }

    /*
     * Redundant full-profile consistency check.
     */
    for (u64 K = 2;
         K <= K_LIMIT;
         ++K) {

        const Witness w =
            profile_winner(r, K);

        if (!w.exists) {
            continue;
        }

        Witness reconstructed;

        for (u64 m = 1; m <= 4; ++m) {
            const Witness candidate =
                best_for_m(
                    r,
                    K,
                    m
                );

            if (!candidate.exists) {
                continue;
            }

            if (
                !reconstructed.exists ||
                candidate.t <
                    reconstructed.t ||
                (
                    candidate.t ==
                    reconstructed.t &&
                    candidate.m <
                    reconstructed.m
                )
            ) {
                reconstructed = candidate;
            }
        }

        if (
            w.m != reconstructed.m ||
            w.k != reconstructed.k ||
            w.t != reconstructed.t
        ) {
            ++stats.profile_failures;
        }
    }
}

static void print_stats(
    const Stats& stats
) {
    std::cout
        << "  CASES="
        << stats.cases
        << '\n';

    std::cout
        << "  WINNER_CHANGES="
        << stats.winner_changes
        << '\n';

    std::cout
        << "  CHANGES_ON_NEW_DIVISOR="
        << stats.changes_on_new_divisor
        << '\n';

    std::cout
        << "  CHANGES_WITHOUT_NEW_DIVISOR="
        << stats.changes_without_new_divisor
        << '\n';

    std::cout
        << "  RECORD_EVENTS="
        << stats.record_events
        << '\n';

    std::cout
        << "  M1_RECORDS="
        << stats.m1_records
        << '\n';

    std::cout
        << "  M2_RECORDS="
        << stats.m2_records
        << '\n';

    std::cout
        << "  M3_RECORDS="
        << stats.m3_records
        << '\n';

    std::cout
        << "  M4_RECORDS="
        << stats.m4_records
        << '\n';

    std::cout
        << "  TRANSITION_1_TO_2="
        << stats.transition_1_to_2
        << '\n';

    std::cout
        << "  TRANSITION_2_TO_3="
        << stats.transition_2_to_3
        << '\n';

    std::cout
        << "  TRANSITION_3_TO_4="
        << stats.transition_3_to_4
        << '\n';

    if (
        stats.min_event_gap !=
        std::numeric_limits<u64>::max()
    ) {
        std::cout
            << "  MIN_EVENT_GAP="
            << stats.min_event_gap
            << '\n';

        std::cout
            << "  MAX_EVENT_GAP="
            << stats.max_event_gap
            << '\n';

        std::cout
            << "  AVERAGE_EVENT_GAP="
            << static_cast<double>(
                   stats.sum_event_gap
               ) /
               static_cast<double>(
                   stats.record_events
               )
            << '\n';
    }

    std::cout
        << "  WINNER_STAYS_BETWEEN_EVENTS="
        << stats.winner_stays_between_events
        << '\n';

    std::cout
        << "  PROFILE_FAILURES="
        << stats.profile_failures
        << '\n';

    std::cout
        << "  VERIFICATION_FAILURES="
        << stats.verification_failures
        << '\n';
}

int main() {
    constexpr u64 EXPERIMENT = 403;
    constexpr u64 CASES = 3000;
    constexpr u64 K_LIMIT = 40;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "K_LIMIT="
        << K_LIMIT
        << '\n';

    std::cout
        << "M_MAX=4\n";

    std::cout
        << "PRIME_MIN=10000\n";

    std::cout
        << "PRIME_MAX=1000000\n\n";

    std::mt19937_64 rng(
        403123456789ULL
    );

    std::vector<PrimePair> cases;
    cases.reserve(CASES);

    for (u64 i = 0;
         i < CASES;
         ++i) {
        cases.push_back(
            generate_case(rng)
        );
    }

    Stats p_stats;
    Stats q_stats;

    for (const PrimePair& c : cases) {
        analyze_prime(
            c.p,
            K_LIMIT,
            p_stats,
            true
        );

        analyze_prime(
            c.q,
            K_LIMIT,
            q_stats,
            false
        );
    }

    std::cout
        << "\nP_STATS\n";

    print_stats(p_stats);

    std::cout
        << "\nQ_STATS\n";

    print_stats(q_stats);

    const bool pass =
        p_stats.profile_failures == 0 &&
        q_stats.profile_failures == 0 &&
        p_stats.verification_failures == 0 &&
        q_stats.verification_failures == 0 &&
        p_stats.changes_without_new_divisor == 0 &&
        q_stats.changes_without_new_divisor == 0;

    std::cout
        << "\nSTATUS="
        << (pass ? "PASS" : "FAIL")
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
