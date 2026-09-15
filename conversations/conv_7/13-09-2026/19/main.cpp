#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
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

std::vector<u64> predicted_starts(
    u64 p,
    const Structure& s
) {
    std::vector<u64> result;

    for (u64 j = 0;
         j < s.q;
         ++j) {

        if (!digitwise_leq(
                p,
                j,
                s.q
            )) {
            continue;
        }

        result.push_back(
            s.s0 +
            j * s.modulus
        );
    }

    return result;
}

std::string vector_string(
    const std::vector<u64>& values,
    std::size_t limit = 30
) {
    std::string s = "[";

    const std::size_t n =
        std::min(
            values.size(),
            limit
        );

    for (std::size_t i = 0;
         i < n;
         ++i) {

        if (i != 0) {
            s += ",";
        }

        s += std::to_string(values[i]);
    }

    if (values.size() > limit) {
        s += ",...";
    }

    s += "]";

    return s;
}

/*
    Check the internal endpoint identity:

        interval[i].hi =
            interval[i+1].lo - 1

    for all non-final intervals.
*/
bool check_internal_endpoints(
    const std::vector<Interval>& intervals
) {
    for (std::size_t i = 0;
         i + 1 < intervals.size();
         ++i) {

        if (intervals[i].hi + 1 !=
            intervals[i + 1].lo) {

            return false;
        }
    }

    return true;
}

/*
    Test whether the final endpoint is m-1,
    m, or something else.
*/
int final_endpoint_type(
    u64 m,
    const std::vector<Interval>& intervals
) {
    if (intervals.empty()) {
        return -1;
    }

    const u64 endpoint =
        intervals.back().hi;

    if (endpoint == m - 1) {
        return 1;
    }

    if (endpoint == m) {
        return 2;
    }

    return 3;
}

void print_failure(
    u64 p,
    u64 m,
    const Structure& s,
    const std::vector<Interval>& actual,
    const std::vector<u64>& starts
) {
    std::cout
        << "\nFAILURE\n"
        << "p=" << p
        << " m=" << m
        << " s0=" << s.s0
        << " modulus=" << s.modulus
        << " q=" << s.q
        << "\n"
        << "predicted_starts="
        << vector_string(starts)
        << "\n"
        << "actual_intervals=";

    for (std::size_t i = 0;
         i < actual.size() &&
         i < 30;
         ++i) {

        if (i != 0) {
            std::cout << " ";
        }

        std::cout
            << "["
            << actual[i].lo
            << ","
            << actual[i].hi
            << "]";
    }

    std::cout << "\n";
}

void run_exhaustive() {
    std::cout
        << "\nPHASE 1: EXACT ENDPOINT EXHAUSTIVE TEST\n";

    u64 total = 0;

    u64 start_set_pass = 0;
    u64 internal_endpoint_pass = 0;

    u64 final_m_minus_1 = 0;
    u64 final_m = 0;
    u64 final_other = 0;

    u64 complete_pass = 0;
    u64 failures = 0;

    for (u64 p = 2;
         p <= 31;
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

            ++total;

            const auto starts =
                predicted_starts(
                    p,
                    s
                );

            bool starts_ok =
                starts.size() ==
                actual.size();

            if (starts_ok) {
                for (std::size_t i = 0;
                     i < actual.size();
                     ++i) {

                    if (actual[i].lo !=
                        starts[i]) {

                        starts_ok = false;
                        break;
                    }
                }
            }

            if (starts_ok) {
                ++start_set_pass;
            }

            const bool internal_ok =
                check_internal_endpoints(
                    actual
                );

            if (internal_ok) {
                ++internal_endpoint_pass;
            }

            const int final_type =
                final_endpoint_type(
                    m,
                    actual
                );

            if (final_type == 1) {
                ++final_m_minus_1;
            } else if (final_type == 2) {
                ++final_m;
            } else if (final_type == 3) {
                ++final_other;
            }

            const bool complete =
                starts_ok &&
                internal_ok &&
                final_type == 1;

            if (complete) {
                ++complete_pass;
            } else {
                ++failures;

                if (failures <= 30) {
                    print_failure(
                        p,
                        m,
                        s,
                        actual,
                        starts
                    );
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 1\n"
        << "multi_interval_cases="
        << total
        << "\n"
        << "start_set_pass="
        << start_set_pass
        << "/" << total
        << "\n"
        << "internal_endpoint_pass="
        << internal_endpoint_pass
        << "/" << total
        << "\n"
        << "final_endpoint_m_minus_1="
        << final_m_minus_1
        << "\n"
        << "final_endpoint_m="
        << final_m
        << "\n"
        << "final_endpoint_other="
        << final_other
        << "\n"
        << "complete_pass="
        << complete_pass
        << "/" << total
        << "\n"
        << "failures="
        << failures
        << "\n";
}

void run_targeted() {
    std::cout
        << "\nPHASE 2: ENDPOINT TARGETS\n";

    struct TestCase {
        u64 p;
        u64 m;
    };

    const std::vector<TestCase> cases = {
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

        if (actual.size() < 2) {
            continue;
        }

        std::cout
            << "p=" << c.p
            << " m=" << c.m
            << " q=" << s.q
            << " intervals=";

        for (const auto& x : actual) {
            std::cout
                << "["
                << x.lo
                << ","
                << x.hi
                << "] ";
        }

        std::cout
            << "\n";
    }
}

void run_large_boundary_sampling() {
    std::cout
        << "\nPHASE 3: LARGE DIRECT BOUNDARY SAMPLING\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7,
        11, 13, 17, 19,
        23, 29, 31, 37,
        43, 53, 67,
        101, 127, 211,
        431, 1009, 2003
    };

    constexpr int CASES = 20000;

    u64 state =
        0x1732026ULL;

    u64 tested = 0;
    u64 pass = 0;
    u64 fail = 0;

    for (int i = 0;
         i < CASES;
         ++i) {

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

        /*
            The internal endpoint relation can be tested
            locally if we have two consecutive valid digitwise
            indices.

            Generate a random valid j and check:

                s(j)+1 ... next_start-1

            for a small sample only.
        */
        if (s.q < 2) {
            continue;
        }

        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        const u64 j =
            state % s.q;

        if (!digitwise_leq(
                p,
                j,
                s.q
            )) {
            continue;
        }

        /*
            Find the next valid subdigit index by
            incrementing only a small number of candidates.
        */
        u64 next_j = j + 1;

        while (next_j < s.q &&
               !digitwise_leq(
                   p,
                   next_j,
                   s.q
               )) {

            ++next_j;
        }

        if (next_j >= s.q) {
            continue;
        }

        const u64 start =
            s.s0 +
            j *
            s.modulus;

        const u64 next_start =
            s.s0 +
            next_j *
            s.modulus;

        if (start > m ||
            next_start > m) {
            continue;
        }

        ++tested;

        /*
            Instead of scanning the interval, test only
            the two crucial endpoints:

                next_start - 1 must be HIT
                next_start must be MISS

            and start must itself be HIT.
        */
        const bool current_hit =
            lucas_hit(
                p,
                m,
                start
            );

        const bool before_next_hit =
            lucas_hit(
                p,
                m,
                next_start - 1
            );

        const bool next_miss =
            !lucas_hit(
                p,
                m,
                next_start
            );

        if (current_hit &&
            before_next_hit &&
            next_miss) {

            ++pass;
        } else {
            ++fail;

            if (fail <= 10) {
                std::cout
                    << "FAILURE "
                    << "p=" << p
                    << " m=" << m
                    << " j=" << j
                    << " next_j=" << next_j
                    << "\n";
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 3\n"
        << "sampled_boundaries="
        << tested
        << "\n"
        << "boundary_pass="
        << pass
        << "/" << tested
        << "\n"
        << "failures="
        << fail
        << "\n";
}

int main() {
    constexpr int EXPERIMENT = 173;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    run_exhaustive();
    run_targeted();
    run_large_boundary_sampling();

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
