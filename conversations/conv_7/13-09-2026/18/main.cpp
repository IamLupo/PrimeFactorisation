#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
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

    for (int i = 0; i < s.a; ++i) {
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

    for (int i = 0; i < s.e; ++i) {
        s.modulus *= p;
    }

    s.q =
        m / s.modulus;

    return s;
}

std::vector<int> base_p_digits(
    u64 x,
    u64 p
) {
    std::vector<int> digits;

    while (x > 0) {
        digits.push_back(
            static_cast<int>(x % p)
        );

        x /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    return digits;
}

std::string digits_string(
    const std::vector<int>& digits
) {
    std::string out = "[";

    for (std::size_t i = 0;
         i < digits.size();
         ++i) {

        if (i != 0) {
            out += ",";
        }

        out += std::to_string(digits[i]);
    }

    out += "]";

    return out;
}

/*
    Construct a random j satisfying

        j <=_p q.

    This avoids enumerating [0,q).
*/
u64 random_subdigit(
    u64 p,
    u64 q,
    u64& state
) {
    std::vector<u64> digits;

    u64 x = q;

    while (x > 0) {
        digits.push_back(x % p);
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

        const u64 chosen =
            state % (digit + 1);

        j += chosen * place;
        place *= p;
    }

    return j;
}

/*
    Construct a random j which is deliberately NOT digitwise
    <= q, whenever possible.
*/
bool random_non_subdigit(
    u64 p,
    u64 q,
    u64& state,
    u64& result
) {
    if (q == 0) {
        return false;
    }

    const auto qdigits =
        base_p_digits(q, p);

    for (int attempt = 0;
         attempt < 32;
         ++attempt) {

        u64 j = 0;
        u64 place = 1;
        bool forced_violation = false;

        for (u64 digit : qdigits) {
            state =
                state *
                6364136223846793005ULL +
                1442695040888963407ULL;

            const u64 random_digit =
                state % p;

            j +=
                random_digit *
                place;

            if (random_digit > digit) {
                forced_violation = true;
            }

            place *= p;
        }

        if (forced_violation &&
            j < q) {

            result = j;
            return true;
        }
    }

    /*
        Deterministic fallback: search the first digit where
        q has room for a violation.
    */
    u64 place = 1;
    u64 x = q;

    while (x > 0) {
        const u64 digit = x % p;

        if (digit + 1 < p) {
            result =
                (digit + 1) * place;

            if (result < q) {
                return true;
            }
        }

        x /= p;
        place *= p;
    }

    return false;
}

/*
    Verify one predicted lattice point as a genuine HIT start.
*/
bool verify_start(
    u64 p,
    u64 m,
    const Structure& s,
    u64 j
) {
    if (j >= s.q) {
        return false;
    }

    if (!digitwise_leq(
            p,
            j,
            s.q
        )) {
        return false;
    }

    const u64 start =
        s.s0 +
        j * s.modulus;

    if (start > m) {
        return false;
    }

    const bool hit =
        lucas_hit(
            p,
            m,
            start
        );

    const bool previous_miss =
        start == 0 ||
        !lucas_hit(
            p,
            m,
            start - 1
        );

    return hit &&
           previous_miss;
}

/*
    Verify one predicted NON-start.
*/
bool verify_non_start(
    u64 p,
    u64 m,
    const Structure& s,
    u64 j
) {
    if (j >= s.q) {
        return false;
    }

    if (digitwise_leq(
            p,
            j,
            s.q
        )) {
        return false;
    }

    const u64 start =
        s.s0 +
        j * s.modulus;

    if (start > m) {
        return false;
    }

    /*
        If the conjecture is correct, this lattice point
        itself should NOT be a HIT start.
    */
    const bool hit =
        lucas_hit(
            p,
            m,
            start
        );

    const bool previous_miss =
        start == 0 ||
        !lucas_hit(
            p,
            m,
            start - 1
        );

    return !(hit && previous_miss);
}

void run_small_exhaustive() {
    std::cout
        << "\nPHASE 1: SMALL COMPLETE VALIDATION\n";

    u64 cases = 0;
    u64 start_tests = 0;
    u64 start_pass = 0;
    u64 nonstart_tests = 0;
    u64 nonstart_pass = 0;

    for (u64 p = 2; p <= 31; ++p) {
        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 1; m <= 1200; ++m) {
            const Structure s =
                analyze_structure(
                    p,
                    m
                );

            if (s.q < 2) {
                continue;
            }

            ++cases;

            /*
                Small q: enumerate all j exactly.
                This is safe because m is small here.
            */
            for (u64 j = 0;
                 j < s.q;
                 ++j) {

                const bool expected =
                    digitwise_leq(
                        p,
                        j,
                        s.q
                    );

                const u64 start =
                    s.s0 +
                    j * s.modulus;

                if (start > m) {
                    continue;
                }

                const bool actual =
                    lucas_hit(
                        p,
                        m,
                        start
                    ) &&
                    !lucas_hit(
                        p,
                        m,
                        start - 1
                    );

                if (expected) {
                    ++start_tests;

                    if (actual) {
                        ++start_pass;
                    } else if (start_pass < 20) {
                        std::cout
                            << "START FAILURE "
                            << "p=" << p
                            << " m=" << m
                            << " q=" << s.q
                            << " j=" << j
                            << "\n";
                    }
                } else {
                    ++nonstart_tests;

                    if (!actual) {
                        ++nonstart_pass;
                    } else if (nonstart_pass < 20) {
                        std::cout
                            << "NONSTART FAILURE "
                            << "p=" << p
                            << " m=" << m
                            << " q=" << s.q
                            << " j=" << j
                            << "\n";
                    }
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 1\n"
        << "cases="
        << cases
        << "\n"
        << "start_predicate_pass="
        << start_pass
        << "/" << start_tests
        << "\n"
        << "nonstart_predicate_pass="
        << nonstart_pass
        << "/" << nonstart_tests
        << "\n";
}

void run_large_random() {
    std::cout
        << "\nPHASE 2: LARGE RANDOM INDEX VALIDATION\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7,
        11, 13, 17, 19,
        23, 29, 31, 37,
        43, 53, 67,
        101, 127, 211,
        431, 1009, 2003
    };

    constexpr int CASES = 10000;
    constexpr int SAMPLES_PER_CASE = 20;

    u64 state =
        0x172B2026ULL;

    u64 start_tests = 0;
    u64 start_pass = 0;

    u64 nonstart_tests = 0;
    u64 nonstart_pass = 0;

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

            /*
                Guaranteed-positive sample.
            */
            const u64 good_j =
                random_subdigit(
                    p,
                    s.q,
                    state
                );

            const u64 good_start =
                s.s0 +
                good_j *
                s.modulus;

            if (good_j < s.q &&
                good_start <= m) {

                ++start_tests;

                if (verify_start(
                        p,
                        m,
                        s,
                        good_j
                    )) {

                    ++start_pass;
                } else if (start_pass < 20) {
                    std::cout
                        << "START FAILURE "
                        << "p=" << p
                        << " m=" << m
                        << " q=" << s.q
                        << " j=" << good_j
                        << "\n";
                }
            }

            /*
                Guaranteed-negative sample when possible.
            */
            u64 bad_j = 0;

            if (random_non_subdigit(
                    p,
                    s.q,
                    state,
                    bad_j
                )) {

                const u64 bad_start =
                    s.s0 +
                    bad_j *
                    s.modulus;

                if (bad_start <= m) {
                    ++nonstart_tests;

                    if (verify_non_start(
                            p,
                            m,
                            s,
                            bad_j
                        )) {

                        ++nonstart_pass;
                    } else if (nonstart_pass < 20) {
                        std::cout
                            << "NONSTART FAILURE "
                            << "p=" << p
                            << " m=" << m
                            << " q=" << s.q
                            << " j=" << bad_j
                            << "\n";
                    }
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 2\n"
        << "start_tests="
        << start_tests
        << "\n"
        << "start_pass="
        << start_pass
        << "/" << start_tests
        << "\n"
        << "nonstart_tests="
        << nonstart_tests
        << "\n"
        << "nonstart_pass="
        << nonstart_pass
        << "/" << nonstart_tests
        << "\n";
}

void run_examples() {
    std::cout
        << "\nPHASE 3: INDEX STRUCTURE EXAMPLES\n";

    struct TestCase {
        u64 p;
        u64 m;
    };

    const std::vector<TestCase> cases = {
        {2, 10},
        {2, 12},
        {2, 20},
        {2, 28},
        {2, 40},
        {3, 41},
        {3, 59},
        {3, 68},
        {3, 125},
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

        std::cout
            << "p=" << c.p
            << " m=" << c.m
            << " e=" << s.e
            << " s0=" << s.s0
            << " p^e=" << s.modulus
            << " q=" << s.q
            << "\n";
    }
}

int main() {
    constexpr int EXPERIMENT = 172;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "B\n";

    run_small_exhaustive();
    run_large_random();
    run_examples();

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "B\n";

    return 0;
}
