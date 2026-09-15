#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using u64 = std::uint64_t;

struct Interval {
    u64 lo;
    u64 hi;
};

struct Structure {
    int a;
    int z;
    int e;

    u64 b;
    u64 s0;
    u64 modulus;
    u64 q;
};

bool is_prime(u64 n) {
    if (n < 2) {
        return false;
    }

    if (n % 2 == 0) {
        return n == 2;
    }

    for (u64 d = 3; d * d <= n; d += 2) {
        if (n % d == 0) {
            return false;
        }
    }

    return true;
}

bool lucas_hit(
    u64 p,
    u64 m,
    u64 t
) {
    while (m > 0 || t > 0) {
        const u64 md = m % p;
        const u64 td = t % p;

        if (td > md) {
            return true;
        }

        m /= p;
        t /= p;
    }

    return false;
}

bool digitwise_leq(
    u64 p,
    u64 a,
    u64 b
) {
    while (a > 0 || b > 0) {
        const u64 ad = a % p;
        const u64 bd = b % p;

        if (ad > bd) {
            return false;
        }

        a /= p;
        b /= p;
    }

    return true;
}

Structure analyze_structure(
    u64 p,
    u64 m
) {
    Structure s{};

    u64 x = m;

    while (x > 0 &&
           x % p == p - 1) {

        ++s.a;
        x /= p;
    }

    s.b = x % p;

    u64 pa = 1;

    for (int i = 0;
         i < s.a;
         ++i) {

        pa *= p;
    }

    x /= p;

    while (x > 0 &&
           x % p == 0) {

        ++s.z;
        x /= p;
    }

    s.e =
        s.a +
        1 +
        s.z;

    s.s0 =
        (s.b + 1) * pa;

    s.modulus = 1;

    for (int i = 0;
         i < s.e;
         ++i) {

        s.modulus *= p;
    }

    s.q =
        m / s.modulus;

    return s;
}

/*
    Next valid digitwise subnumber.

    Input:
        j <=_p q

    Output:
        smallest k > j such that k <=_p q.

    The algorithm searches from the least-significant digit
    upward for the first place where j_i < q_i.

    That digit can be increased by one, and all lower digits
    are set to zero.
*/
bool next_valid_index(
    u64 p,
    u64 j,
    u64 q,
    u64& result
) {
    u64 xj = j;
    u64 xq = q;

    u64 place = 1;

    while (xj > 0 ||
           xq > 0) {

        const u64 jd = xj % p;
        const u64 qd = xq % p;

        /*
            This position can be increased while remaining
            digitwise <= q.
        */
        if (jd < qd) {
            /*
                Keep all higher digits of j unchanged.
            */
            const u64 lower_mask =
                place - 1;

            const u64 higher_part =
                j -
                (j % (place * p));

            result =
                higher_part +
                (jd + 1) * place;

            /*
                All lower digits become zero.
            */
            (void)lower_mask;

            return true;
        }

        xj /= p;
        xq /= p;
        place *= p;
    }

    return false;
}

/*
    Generate the complete set of valid indices for small q.
*/
std::vector<u64> all_valid_indices(
    u64 p,
    u64 q
) {
    std::vector<u64> result;

    for (u64 j = 0;
         j < q;
         ++j) {

        if (digitwise_leq(
                p,
                j,
                q
            )) {

            result.push_back(j);
        }
    }

    return result;
}

/*
    Actual intervals for exhaustive validation only.
*/
std::vector<Interval> actual_intervals(
    u64 p,
    u64 m
) {
    std::vector<Interval> result;

    bool inside = false;
    u64 start = 0;

    for (u64 t = 1;
         t <= m;
         ++t) {

        const bool hit =
            lucas_hit(
                p,
                m,
                t
            );

        if (hit && !inside) {
            start = t;
            inside = true;
        }

        if (!hit && inside) {
            result.push_back({
                start,
                t - 1
            });

            inside = false;
        }
    }

    if (inside) {
        result.push_back({
            start,
            m
        });
    }

    return result;
}

/*
    Construct the predicted interval for a valid index j.

    j_next is allowed to equal q. In that case

        endpoint = q*p^e - 1.

    Since

        m = q*p^e + s0 - 1,

    the final MISS region follows immediately after the
    predicted final interval.
*/
Interval predicted_interval(
    u64 p,
    const Structure& s,
    u64 j
) {
    u64 next = 0;

    const bool found =
        next_valid_index(
            p,
            j,
            s.q,
            next
        );

    /*
        For every j < q which is valid, q itself is always
        digitwise <= q, so a next index must exist.
    */
    if (!found) {
        next = s.q;
    }

    return {
        s.s0 +
        j * s.modulus,

        next * s.modulus -
        1
    };
}

std::string intervals_string(
    const std::vector<Interval>& intervals,
    std::size_t limit = 30
) {
    std::string result;

    const std::size_t n =
        std::min(
            intervals.size(),
            limit
        );

    for (std::size_t i = 0;
         i < n;
         ++i) {

        if (i != 0) {
            result += " ";
        }

        result +=
            "[" +
            std::to_string(
                intervals[i].lo
            ) +
            "," +
            std::to_string(
                intervals[i].hi
            ) +
            "]";
    }

    if (intervals.size() > limit) {
        result += " ...";
    }

    return result;
}

void run_exhaustive_phase() {
    std::cout
        << "\nPHASE 1: COMPLETE CLOSED-FORM INTERVAL TEST\n";

    u64 total_cases = 0;
    u64 total_intervals = 0;

    u64 index_transition_pass = 0;
    u64 interval_formula_pass = 0;

    u64 failures = 0;

    for (u64 p = 2;
         p <= 37;
         ++p) {

        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 1;
             m <= 2000;
             ++m) {

            const Structure s =
                analyze_structure(
                    p,
                    m
                );

            const auto actual =
                actual_intervals(
                    p,
                    m
                );

            if (actual.size() < 2) {
                continue;
            }

            ++total_cases;
            total_intervals +=
                actual.size();

            const auto valid =
                all_valid_indices(
                    p,
                    s.q
                );

            bool case_ok = true;

            /*
                Verify every interval independently.
            */
            if (valid.size() !=
                actual.size()) {

                case_ok = false;
            }

            const std::size_t count =
                std::min(
                    valid.size(),
                    actual.size()
                );

            for (std::size_t i = 0;
                 i < count;
                 ++i) {

                const Interval predicted =
                    predicted_interval(
                        p,
                        s,
                        valid[i]
                    );

                if (predicted.lo !=
                        actual[i].lo ||
                    predicted.hi !=
                        actual[i].hi) {

                    case_ok = false;

                    if (failures < 20) {
                        std::cout
                            << "\nFAILURE\n"
                            << "p=" << p
                            << " m=" << m
                            << " j=" << valid[i]
                            << "\n"
                            << "predicted=["
                            << predicted.lo
                            << ","
                            << predicted.hi
                            << "]\n"
                            << "actual=["
                            << actual[i].lo
                            << ","
                            << actual[i].hi
                            << "]\n"
                            << "actual_all="
                            << intervals_string(
                                actual
                            )
                            << "\n";
                    }

                    ++failures;
                }
            }

            /*
                Verify the index transition independently.
            */
            bool transition_ok = true;

            for (std::size_t i = 0;
                 i + 1 < valid.size();
                 ++i) {

                u64 next = 0;

                if (!next_valid_index(
                        p,
                        valid[i],
                        s.q,
                        next
                    )) {

                    transition_ok = false;
                    break;
                }

                if (next !=
                    valid[i + 1]) {

                    transition_ok = false;
                    break;
                }
            }

            if (transition_ok) {
                ++index_transition_pass;
            }

            if (case_ok) {
                ++interval_formula_pass;
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 1\n"
        << "multi_interval_cases="
        << total_cases
        << "\n"
        << "actual_intervals="
        << total_intervals
        << "\n"
        << "index_transition_pass="
        << index_transition_pass
        << "/"
        << total_cases
        << "\n"
        << "complete_interval_formula_pass="
        << interval_formula_pass
        << "/"
        << total_cases
        << "\n"
        << "failures="
        << failures
        << "\n";
}

void run_targeted_phase() {
    std::cout
        << "\nPHASE 2: TARGETED CLOSED-FORM EXAMPLES\n";

    struct TestCase {
        u64 p;
        u64 m;
    };

    const std::vector<TestCase> cases = {
        {2, 6},
        {2, 10},
        {2, 12},
        {2, 20},
        {2, 24},
        {2, 28},
        {2, 40},
        {2, 48},

        {3, 22},
        {3, 41},
        {3, 59},
        {3, 68},
        {3, 125},
        {3, 377},

        {5, 59},
        {5, 94},
        {5, 194},
        {5, 219}
    };

    for (const auto& c : cases) {
        const Structure s =
            analyze_structure(
                c.p,
                c.m
            );

        const auto actual =
            actual_intervals(
                c.p,
                c.m
            );

        const auto valid =
            all_valid_indices(
                c.p,
                s.q
            );

        std::cout
            << "\np="
            << c.p
            << " m="
            << c.m
            << "\n";

        std::cout
            << "q="
            << s.q
            << " s0="
            << s.s0
            << " p^e="
            << s.modulus
            << "\n";

        std::cout
            << "actual="
            << intervals_string(
                actual
            )
            << "\n";

        std::cout
            << "formula=";

        for (u64 j : valid) {
            const Interval interval =
                predicted_interval(
                    c.p,
                    s,
                    j
                );

            std::cout
                << "["
                << interval.lo
                << ","
                << interval.hi
                << "] ";
        }

        std::cout
            << "\n";
    }
}

/*
    Large validation.

    We never enumerate q or the complete index set.

    Instead we randomly construct a valid j and calculate
    its successor directly.
*/
u64 random_valid_index(
    u64 p,
    u64 q,
    u64& state
) {
    std::vector<u64> digits;

    u64 x = q;

    while (x > 0) {
        digits.push_back(
            x % p
        );

        x /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    u64 j = 0;
    u64 place = 1;

    for (u64 digit : digits) {
        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        const u64 selected =
            state % (digit + 1);

        j +=
            selected *
            place;

        place *= p;
    }

    return j;
}

void run_large_phase() {
    std::cout
        << "\nPHASE 3: LARGE DIRECT FORMULA VALIDATION\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7,
        11, 13, 17, 19,
        23, 29, 31, 37,
        43, 53, 67,
        101, 127, 211,
        431, 1009,
        2003, 4001
    };

    constexpr int CASES = 50000;
    constexpr int SAMPLES_PER_CASE = 10;

    u64 state =
        0x1752026ULL;

    u64 tested = 0;
    u64 pass = 0;
    u64 fail = 0;

    for (int c = 0;
         c < CASES;
         ++c) {

        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        const u64 p =
            primes[
                state %
                primes.size()
            ];

        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        const u64 m =
            1 +
            state %
            1000000000000000000ULL;

        const Structure s =
            analyze_structure(
                p,
                m
            );

        if (s.q < 2) {
            continue;
        }

        for (int k = 0;
             k < SAMPLES_PER_CASE;
             ++k) {

            const u64 j =
                random_valid_index(
                    p,
                    s.q,
                    state
                );

            const u64 start =
                s.s0 +
                j *
                s.modulus;

            if (start > m) {
                continue;
            }

            u64 next = 0;

            if (!next_valid_index(
                    p,
                    j,
                    s.q,
                    next
                )) {

                continue;
            }

            const Interval predicted =
                predicted_interval(
                    p,
                    s,
                    j
                );

            /*
                Direct boundary validation:

                    start     = HIT
                    start - 1 = MISS
                    endpoint  = HIT
                    endpoint+1 = MISS

                No scan over the interval.
            */
            const bool start_ok =
                lucas_hit(
                    p,
                    m,
                    predicted.lo
                ) &&
                predicted.lo > 0 &&
                !lucas_hit(
                    p,
                    m,
                    predicted.lo - 1
                );

            const bool end_ok =
                lucas_hit(
                    p,
                    m,
                    predicted.hi
                ) &&
                (
                    predicted.hi >= m ||
                    !lucas_hit(
                        p,
                        m,
                        predicted.hi + 1
                    )
                );

            ++tested;

            if (start_ok &&
                end_ok) {

                ++pass;
            } else {
                ++fail;

                if (fail <= 20) {
                    std::cout
                        << "FAILURE "
                        << "p=" << p
                        << " m=" << m
                        << " j=" << j
                        << " next="
                        << next
                        << " interval=["
                        << predicted.lo
                        << ","
                        << predicted.hi
                        << "]\n";
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 3\n"
        << "sampled_intervals="
        << tested
        << "\n"
        << "boundary_formula_pass="
        << pass
        << "/"
        << tested
        << "\n"
        << "failures="
        << fail
        << "\n";
}

int main() {
    constexpr int EXPERIMENT = 175;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    run_exhaustive_phase();
    run_targeted_phase();
    run_large_phase();

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
