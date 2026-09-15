#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Parameters {
    u64 p;
    u64 m;
    u64 a;
    u64 b;
    u64 z;
    u64 e;
    u64 s0;
    u64 q;
    u64 pe;
};

struct Interval {
    u64 start;
    u64 end;
};

void print_start() {
    std::cout << "START EXPERIMENT 177\n";
}

void print_finish() {
    std::cout << "FINISHED EXPERIMENT 177\n";
}

u64 pow_u64(u64 base, u64 exp) {
    u128 r = 1;

    for (u64 i = 0; i < exp; ++i) {
        r *= base;
    }

    return static_cast<u64>(r);
}

bool digitwise_leq(u64 p, u64 j, u64 q) {
    while (j > 0 || q > 0) {
        u64 jd = j % p;
        u64 qd = q % p;

        if (jd > qd) {
            return false;
        }

        j /= p;
        q /= p;
    }

    return true;
}

bool is_hit(u64 p, u64 m, u64 t) {
    u64 a = t;
    u64 b = m;

    while (a > 0 || b > 0) {
        u64 ad = a % p;
        u64 bd = b % p;

        if (ad > bd) {
            return true;
        }

        a /= p;
        b /= p;
    }

    return false;
}

bool derive_parameters(u64 p, u64 m, Parameters &out) {
    u64 x = m;
    u64 a = 0;

    while (x > 0 && x % p == p - 1) {
        ++a;
        x /= p;
    }

    if (x == 0) {
        return false;
    }

    u64 b = x % p;
    x /= p;

    if (b >= p - 1) {
        return false;
    }

    u64 z = 0;

    while (x > 0 && x % p == 0) {
        ++z;
        x /= p;
    }

    u64 e = a + 1 + z;
    u64 pa = pow_u64(p, a);
    u64 pe = pow_u64(p, e);

    u128 numerator =
        static_cast<u128>(m) + 1 -
        static_cast<u128>((b + 1) * pa);

    if (numerator == 0 || numerator % pe != 0) {
        return false;
    }

    u64 q = static_cast<u64>(numerator / pe);

    if (q == 0) {
        return false;
    }

    out = {
        p,
        m,
        a,
        b,
        z,
        e,
        (b + 1) * pa,
        q,
        pe
    };

    return true;
}

u64 lowest_available_digit(u64 p, u64 j, u64 q) {
    u64 k = 0;

    while (j > 0 || q > 0) {
        u64 jd = j % p;
        u64 qd = q % p;

        if (jd < qd) {
            return k;
        }

        j /= p;
        q /= p;
        ++k;
    }

    return UINT64_MAX;
}

u64 next_valid_index(u64 p, u64 j, u64 q) {
    u64 k = lowest_available_digit(p, j, q);

    if (k == UINT64_MAX) {
        return q;
    }

    u64 pk = pow_u64(p, k);

    return j + pk - (j % pk);
}

u64 predicted_gap(u64 p, u64 j, u64 q) {
    u64 jp = next_valid_index(p, j, q);
    return jp - j;
}

std::vector<Interval> actual_intervals(u64 p, u64 m) {
    std::vector<Interval> result;

    bool inside = false;
    u64 start = 0;

    for (u64 t = 1; t <= m; ++t) {
        bool hit = is_hit(p, m, t);

        if (hit && !inside) {
            inside = true;
            start = t;
        }

        if (!hit && inside) {
            inside = false;
            result.push_back({start, t - 1});
        }
    }

    if (inside) {
        result.push_back({start, m});
    }

    return result;
}

std::vector<Interval> predicted_intervals(const Parameters &par) {
    std::vector<Interval> result;

    for (u64 j = 0; j < par.q; ++j) {
        if (!digitwise_leq(par.p, j, par.q)) {
            continue;
        }

        u64 jp = next_valid_index(par.p, j, par.q);

        u128 start =
            static_cast<u128>(par.s0) +
            static_cast<u128>(j) * par.pe;

        u128 end =
            static_cast<u128>(par.s0) +
            static_cast<u128>(jp) * par.pe -
            2;

        result.push_back({
            static_cast<u64>(start),
            static_cast<u64>(end)
        });
    }

    return result;
}

bool same_intervals(
    const std::vector<Interval> &a,
    const std::vector<Interval> &b
) {
    if (a.size() != b.size()) {
        return false;
    }

    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].start != b[i].start ||
            a[i].end != b[i].end) {
            return false;
        }
    }

    return true;
}

u64 random_valid_index(
    u64 p,
    u64 q,
    std::mt19937_64 &rng
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

    for (;;) {
        u128 value = 0;
        u128 place = 1;

        for (u64 d : digits) {
            std::uniform_int_distribution<u64> dist(0, d);
            u64 chosen = dist(rng);

            value += static_cast<u128>(chosen) * place;
            place *= p;
        }

        u64 j = static_cast<u64>(value);

        if (j < q) {
            return j;
        }
    }
}

u64 random_large_m(std::mt19937_64 &rng) {
    std::uniform_int_distribution<u64> dist(
        1000000000000ULL,
        1000000000000000000ULL
    );

    return dist(rng);
}

void run_small_test(
    u64 &cases,
    u64 &interval_cases,
    u64 &formula_pass,
    u64 &length_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    for (u64 p : primes) {
        for (u64 m = 2; m <= 2000; ++m) {
            Parameters par{};

            if (!derive_parameters(p, m, par)) {
                continue;
            }

            ++cases;

            auto actual = actual_intervals(p, m);
            auto predicted = predicted_intervals(par);

            if (same_intervals(actual, predicted)) {
                ++formula_pass;
            }

            for (u64 j = 0; j < par.q; ++j) {
                if (!digitwise_leq(p, j, par.q)) {
                    continue;
                }

                ++interval_cases;

                u64 jp = next_valid_index(p, j, par.q);
                u64 gap = jp - j;

                u128 expected_length =
                    static_cast<u128>(gap) *
                    par.pe - 1;

                u64 predicted_length =
                    static_cast<u64>(expected_length);

                u64 predicted_start =
                    static_cast<u64>(
                        static_cast<u128>(par.s0) +
                        static_cast<u128>(j) * par.pe
                    );

                u64 predicted_end =
                    predicted_start +
                    predicted_length - 1;

                bool ok =
                    is_hit(p, m, predicted_start) &&
                    is_hit(p, m, predicted_end);

                if (predicted_start > 1) {
                    if (is_hit(
                        p,
                        m,
                        predicted_start - 1
                    )) {
                        ok = false;
                    }
                }

                if (predicted_end < m) {
                    if (is_hit(
                        p,
                        m,
                        predicted_end + 1
                    )) {
                        ok = false;
                    }
                }

                if (ok) {
                    ++length_pass;
                }
            }
        }
    }
}

void run_large_test(
    u64 trials,
    u64 &formula_pass,
    u64 &boundary_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    std::mt19937_64 rng(0x17720260913ULL);

    for (u64 i = 0; i < trials; ++i) {
        u64 p = primes[i % primes.size()];
        u64 m = random_large_m(rng);

        Parameters par{};

        if (!derive_parameters(p, m, par)) {
            --i;
            continue;
        }

        u64 j = random_valid_index(
            p,
            par.q,
            rng
        );

        u64 jp = next_valid_index(
            p,
            j,
            par.q
        );

        u64 gap = jp - j;

        u128 start128 =
            static_cast<u128>(par.s0) +
            static_cast<u128>(j) * par.pe;

        u128 end128 =
            static_cast<u128>(par.s0) +
            static_cast<u128>(jp) * par.pe -
            2;

        u64 start = static_cast<u64>(start128);
        u64 end = static_cast<u64>(end128);

        u64 expected_length =
            static_cast<u64>(
                static_cast<u128>(gap) * par.pe - 1
            );

        if (end - start + 1 == expected_length) {
            ++formula_pass;
        }

        bool ok =
            start <= m &&
            end <= m &&
            start <= end &&
            is_hit(p, m, start) &&
            is_hit(p, m, end);

        if (start > 1 &&
            is_hit(p, m, start - 1)) {
            ok = false;
        }

        if (end < m &&
            is_hit(p, m, end + 1)) {
            ok = false;
        }

        if (ok) {
            ++boundary_pass;
        }
    }
}

int main() {
    print_start();

    u64 cases = 0;
    u64 interval_cases = 0;
    u64 formula_pass = 0;
    u64 length_pass = 0;

    run_small_test(
        cases,
        interval_cases,
        formula_pass,
        length_pass
    );

    std::cout << "\nSMALL EXHAUSTIVE\n";

    std::cout
        << "parameter_cases="
        << cases
        << "\n";

    std::cout
        << "interval_formula_pass="
        << formula_pass
        << "/" << cases
        << "\n";

    std::cout
        << "interval_boundary_cases="
        << interval_cases
        << "\n";

    std::cout
        << "length_boundary_pass="
        << length_pass
        << "/" << interval_cases
        << "\n";

    u64 large_formula_pass = 0;
    u64 large_boundary_pass = 0;

    run_large_test(
        200000,
        large_formula_pass,
        large_boundary_pass
    );

    std::cout << "\nLARGE RANDOM\n";

    std::cout
        << "length_formula_pass="
        << large_formula_pass
        << "/200000\n";

    std::cout
        << "boundary_pass="
        << large_boundary_pass
        << "/200000\n";

    std::cout << "\nCORRECTED IDENTITY\n";

    std::cout
        << "j_next - j = p^k - (j mod p^k)\n";

    std::cout
        << "|I_j| = (j_next - j) * p^e - 1\n";

    std::cout
        << "I_j = [s0 + j*p^e, "
        << "s0 + j_next*p^e - 2]\n";

    print_finish();

    return 0;
}
