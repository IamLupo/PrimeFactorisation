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
    u64 x,
    u64 m
) {
    while (x > 0 || m > 0) {
        const u64 xd = x % p;
        const u64 md = m % p;

        if (xd > md) {
            return false;
        }

        x /= p;
        m /= p;
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
    Actual intervals.

    Used only for small exhaustive validation.
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
    Find the smallest MISS y > x.

    A MISS is a number whose base-p digits are all <=
    the corresponding digits of m.

    Algorithm:

      1. Find the highest digit where x_i > m_i.
      2. At or below that point, x cannot be repaired by
         changing a lower digit.
      3. Find the first higher position k where x_k < m_k.
      4. Increase x_k by one and set all lower digits to 0.
      5. Keep all higher digits unchanged.

    If no such k exists, there is no later MISS.
*/
bool next_miss(
    u64 p,
    u64 m,
    u64 x,
    u64& result
) {
    std::vector<u64> xd;
    std::vector<u64> md;

    u64 tx = x;
    u64 tm = m;

    while (tx > 0 ||
           tm > 0) {

        xd.push_back(tx % p);
        md.push_back(tm % p);

        tx /= p;
        tm /= p;
    }

    while (xd.size() < md.size()) {
        xd.push_back(0);
    }

    while (md.size() < xd.size()) {
        md.push_back(0);
    }

    const int digits =
        static_cast<int>(xd.size());

    int highest_violation = -1;

    for (int i = 0;
         i < digits;
         ++i) {

        if (xd[i] > md[i]) {
            highest_violation = i;
        }
    }

    /*
        If x is already a MISS, the successor is the ordinary
        next digitwise-subnumber. The same construction below
        handles it by taking the first incrementable position.
    */

    int start_position =
        std::max(
            0,
            highest_violation
        );

    /*
        We need a position strictly above every violating
        digit. If there is no violation, position 0 is allowed.
    */
    if (highest_violation >= 0) {
        start_position =
            highest_violation + 1;
    }

    for (int k = start_position;
         k < digits;
         ++k) {

        /*
            All higher digits must already be valid.
        */
        bool higher_valid = true;

        for (int i = k + 1;
             i < digits;
             ++i) {

            if (xd[i] > md[i]) {
                higher_valid = false;
                break;
            }
        }

        if (!higher_valid) {
            continue;
        }

        /*
            We need to increase this digit.
        */
        if (xd[k] >= md[k]) {
            continue;
        }

        std::vector<u64> yd =
            xd;

        yd[k] =
            xd[k] + 1;

        for (int i = 0;
             i < k;
             ++i) {

            yd[i] = 0;
        }

        /*
            Reconstruct.
        */
        u64 value = 0;
        u64 place = 1;

        for (int i = 0;
             i < digits;
             ++i) {

            value +=
                yd[i] * place;

            place *= p;
        }

        if (value > x &&
            digitwise_leq(
                p,
                value,
                m
            )) {

            result = value;
            return true;
        }
    }

    return false;
}

/*
    Predict the endpoint of a HIT interval beginning at start.

    Endpoint is:

        next_miss(start) - 1

    or m if no later MISS exists.
*/
u64 predicted_endpoint(
    u64 p,
    u64 m,
    u64 start
) {
    u64 next = 0;

    if (next_miss(
            p,
            m,
            start,
            next
        )) {

        return next - 1;
    }

    return m;
}

std::string interval_string(
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

void run_small_exhaustive() {
    std::cout
        << "\nPHASE 1: EXACT NEXT-MISS EXHAUSTIVE TEST\n";

    u64 total_intervals = 0;
    u64 endpoint_pass = 0;
    u64 endpoint_fail = 0;

    for (u64 p = 2;
         p <= 37;
         ++p) {

        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 1;
             m <= 2000;
             ++m) {

            const auto actual =
                actual_intervals(
                    p,
                    m
                );

            if (actual.empty()) {
                continue;
            }

            for (const auto& interval :
                 actual) {

                ++total_intervals;

                const u64 predicted =
                    predicted_endpoint(
                        p,
                        m,
                        interval.lo
                    );

                if (predicted ==
                    interval.hi) {

                    ++endpoint_pass;
                } else {
                    ++endpoint_fail;

                    if (endpoint_fail <= 30) {
                        std::cout
                            << "\nFAILURE\n"
                            << "p=" << p
                            << " m=" << m
                            << "\n"
                            << "start="
                            << interval.lo
                            << "\n"
                            << "actual_endpoint="
                            << interval.hi
                            << "\n"
                            << "predicted_endpoint="
                            << predicted
                            << "\n";
                    }
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 1\n"
        << "intervals_tested="
        << total_intervals
        << "\n"
        << "endpoint_pass="
        << endpoint_pass
        << "/"
        << total_intervals
        << "\n"
        << "endpoint_fail="
        << endpoint_fail
        << "\n";
}

void run_targeted() {
    std::cout
        << "\nPHASE 2: TARGETED INTERVALS\n";

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
        const auto actual =
            actual_intervals(
                c.p,
                c.m
            );

        std::cout
            << "\np="
            << c.p
            << " m="
            << c.m
            << "\n";

        std::cout
            << "actual="
            << interval_string(actual)
            << "\n";

        std::cout
            << "predicted=";

        for (const auto& interval :
             actual) {

            const u64 endpoint =
                predicted_endpoint(
                    c.p,
                    c.m,
                    interval.lo
                );

            std::cout
                << "["
                << interval.lo
                << ","
                << endpoint
                << "] ";
        }

        std::cout
            << "\n";
    }
}

/*
    Large direct test.

    We generate random m up to 1e18 and random t that is known
    to be a HIT-start by constructing j <=_p q.

    Then we calculate the predicted endpoint from digits and
    verify only the boundary:

        endpoint     = HIT
        endpoint + 1 = MISS

    No scan over the interval is performed.
*/
bool random_subdigit(
    u64 p,
    u64 q,
    u64& state,
    u64& result
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

    result = 0;

    u64 place = 1;

    for (u64 digit : digits) {
        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        const u64 chosen =
            state %
            (digit + 1);

        result +=
            chosen *
            place;

        place *= p;
    }

    return true;
}

void run_large_direct() {
    std::cout
        << "\nPHASE 3: LARGE DIRECT ENDPOINT TEST\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7,
        11, 13, 17, 19,
        23, 29, 31, 37,
        43, 53, 67,
        101, 127, 211,
        431, 1009, 2003
    };

    constexpr int CASES = 50000;

    u64 state =
        0x1742026ULL;

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

        if (s.q == 0) {
            continue;
        }

        u64 j = 0;

        random_subdigit(
            p,
            s.q,
            state,
            j
        );

        const u64 start =
            s.s0 +
            j *
            s.modulus;

        if (start > m ||
            start == 0) {

            continue;
        }

        const u64 endpoint =
            predicted_endpoint(
                p,
                m,
                start
            );

        /*
            Check endpoint directly.
        */
        const bool endpoint_hit =
            lucas_hit(
                p,
                m,
                endpoint
            );

        const bool after_miss =
            endpoint >= m ||
            !lucas_hit(
                p,
                m,
                endpoint + 1
            );

        const bool start_hit =
            lucas_hit(
                p,
                m,
                start
            );

        const bool start_previous_miss =
            start == 0 ||
            !lucas_hit(
                p,
                m,
                start - 1
            );

        ++tested;

        const bool ok =
            start_hit &&
            start_previous_miss &&
            endpoint_hit &&
            after_miss;

        if (ok) {
            ++pass;
        } else {
            ++fail;

            if (fail <= 20) {
                std::cout
                    << "FAILURE "
                    << "p=" << p
                    << " m=" << m
                    << " start="
                    << start
                    << " endpoint="
                    << endpoint
                    << "\n";
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 3\n"
        << "tested="
        << tested
        << "\n"
        << "endpoint_boundary_pass="
        << pass
        << "/"
        << tested
        << "\n"
        << "failures="
        << fail
        << "\n";
}

int main() {
    constexpr int EXPERIMENT = 174;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    run_small_exhaustive();
    run_targeted();
    run_large_direct();

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
