#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

struct CaseData {
    u64 p;
    u64 q;
};

struct Polynomial {
    std::vector<u64> coeff;
};

struct RootWitness {
    bool found = false;
    u64 root = 0;
    u64 gcd_coeff = 0;
};

struct Factor {
    u64 prime;
    int exponent;
};

std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }

    bool negative = (x < 0);

    u128 v =
        negative
            ? static_cast<u128>(-x)
            : static_cast<u128>(x);

    std::string s;

    while (v > 0) {
        int digit =
            static_cast<int>(v % 10);

        s.push_back(
            static_cast<char>('0' + digit)
        );

        v /= 10;
    }

    if (negative) {
        s.push_back('-');
    }

    std::reverse(s.begin(), s.end());

    return s;
}

u64 gcd_u64(
    u64 a,
    u64 b
) {
    return std::gcd(a, b);
}

u64 integer_sqrt(
    u64 n
) {
    if (n == 0) {
        return 0;
    }

    u64 lo = 0;
    u64 hi =
        std::min<u64>(
            n,
            1ULL << 32
        );

    while (lo + 1 < hi) {
        u64 mid =
            lo + (hi - lo) / 2;

        if (mid <= n / mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    if (hi <= n / hi) {
        return hi;
    }

    return lo;
}

std::vector<u64> base_digits(
    u64 x,
    u64 base
) {
    std::vector<u64> d;

    if (x == 0) {
        d.push_back(0);
        return d;
    }

    while (x > 0) {
        d.push_back(
            x % base
        );

        x /= base;
    }

    return d;
}

u64 miss_prefix(
    u64 m,
    i128 y_signed,
    u64 base
) {
    if (y_signed < 0) {
        return 0;
    }

    u64 y =
        static_cast<u64>(
            y_signed
        );

    std::vector<u64> md =
        base_digits(
            m,
            base
        );

    std::vector<u64> yd =
        base_digits(
            y,
            base
        );

    std::size_t L =
        std::max(
            md.size(),
            yd.size()
        );

    md.resize(L, 0);
    yd.resize(L, 0);

    std::vector<u128> suffix(
        L + 1,
        1
    );

    for (
        std::size_t i = L;
        i-- > 0;
    ) {
        suffix[i] =
            suffix[i + 1] *
            static_cast<u128>(
                md[i] + 1
            );
    }

    u128 answer = 0;
    bool tight = true;

    for (
        std::size_t pos = L;
        pos-- > 0;
    ) {
        if (!tight) {
            break;
        }

        u64 yi = yd[pos];
        u64 mi = md[pos];

        u64 upper =
            std::min(
                yi,
                mi + 1
            );

        for (
            u64 xi = 0;
            xi < upper;
            ++xi
        ) {
            answer += suffix[pos + 1];
        }

        if (yi <= mi) {
            tight = true;
        } else {
            tight = false;
            break;
        }
    }

    if (tight) {
        ++answer;
    }

    return static_cast<u64>(
        answer
    );
}

u64 local_delta(
    u64 old_m,
    u64 new_m,
    i128 y,
    u64 base
) {
    u64 a =
        miss_prefix(
            old_m,
            y,
            base
        );

    u64 b =
        miss_prefix(
            new_m,
            y,
            base
        );

    if (b < a) {
        std::cerr
            << "ERROR: negative derivative\n";

        std::exit(1);
    }

    return b - a;
}

u64 base_power(
    u64 base,
    std::size_t exponent
) {
    u64 result = 1;

    for (
        std::size_t i = 0;
        i < exponent;
        ++i
    ) {
        result *= base;
    }

    return result;
}

Polynomial build_polynomial(
    u64 target_m,
    i128 y,
    u64 base
) {
    std::vector<u64> digits =
        base_digits(
            target_m,
            base
        );

    Polynomial P;

    P.coeff.assign(
        digits.size(),
        0
    );

    u64 current = 0;

    for (
        std::size_t r = 0;
        r < digits.size();
        ++r
    ) {
        u64 step =
            base_power(
                base,
                r
            );

        for (
            u64 k = 0;
            k < digits[r];
            ++k
        ) {
            u64 next =
                current + step;

            u64 delta =
                local_delta(
                    current,
                    next,
                    y,
                    base
                );

            P.coeff[r] += delta;

            current = next;
        }
    }

    if (current != target_m) {
        std::cerr
            << "ERROR: path mismatch\n";

        std::exit(1);
    }

    return P;
}

u64 coefficient_gcd(
    const Polynomial &P
) {
    u64 g = 0;

    for (u64 c : P.coeff) {
        g = gcd_u64(
            g,
            c
        );
    }

    return g;
}

std::vector<Factor> factor_integer(
    u64 n
) {
    std::vector<Factor> factors;

    if (n <= 1) {
        return factors;
    }

    int exponent = 0;

    while ((n & 1ULL) == 0) {
        n >>= 1ULL;
        ++exponent;
    }

    if (exponent != 0) {
        factors.push_back({
            2,
            exponent
        });
    }

    for (
        u64 p = 3;
        p <= n / p;
        p += 2
    ) {
        if (n % p != 0) {
            continue;
        }

        exponent = 0;

        while (n % p == 0) {
            n /= p;
            ++exponent;
        }

        factors.push_back({
            p,
            exponent
        });
    }

    if (n > 1) {
        factors.push_back({
            n,
            1
        });
    }

    return factors;
}

std::string factor_string(
    u64 n
) {
    if (n == 0) {
        return "0";
    }

    if (n == 1) {
        return "1";
    }

    std::vector<Factor> factors =
        factor_integer(n);

    std::string s;

    for (
        std::size_t i = 0;
        i < factors.size();
        ++i
    ) {
        if (i != 0) {
            s += " * ";
        }

        s +=
            std::to_string(
                factors[i].prime
            );

        if (factors[i].exponent > 1) {
            s += "^";
            s +=
                std::to_string(
                    factors[i].exponent
                );
        }
    }

    return s;
}

bool all_coeff_divisible(
    const Polynomial &P,
    u64 d
) {
    if (d == 0) {
        return false;
    }

    for (u64 c : P.coeff) {
        if (c % d != 0) {
            return false;
        }
    }

    return true;
}

u64 mod_add(
    u64 a,
    u64 b,
    u64 mod
) {
    u64 r = a + b;

    if (r >= mod || r < a) {
        r %= mod;
    }

    return r;
}

u64 mod_mul(
    u64 a,
    u64 b,
    u64 mod
) {
    return static_cast<u64>(
        (
            static_cast<u128>(a) *
            static_cast<u128>(b)
        ) %
        static_cast<u128>(mod)
    );
}

u64 evaluate_mod(
    const Polynomial &P,
    u64 x,
    u64 mod
) {
    u64 result = 0;

    for (
        std::size_t i = P.coeff.size();
        i-- > 0;
    ) {
        result =
            mod_mul(
                result,
                x,
                mod
            );

        result =
            mod_add(
                result,
                P.coeff[i] % mod,
                mod
            );
    }

    return result;
}

RootWitness find_first_root(
    const Polynomial &P,
    u64 mod,
    u64 gcd_coeff
) {
    RootWitness W;

    if (gcd_coeff != 1) {
        for (u64 x = 0; x < mod; ++x) {
            if (
                evaluate_mod(
                    P,
                    x,
                    mod
                ) == 0
            ) {
                W.found = true;
                W.root = x;
                W.gcd_coeff = gcd_coeff;
                return W;
            }
        }
    } else {
        for (u64 x = 0; x < mod; ++x) {
            if (
                evaluate_mod(
                    P,
                    x,
                    mod
                ) == 0
            ) {
                W.found = true;
                W.root = x;
                W.gcd_coeff = 1;
                return W;
            }
        }
    }

    return W;
}

void print_coefficients(
    const Polynomial &P,
    const std::string &prefix
) {
    std::cout
        << prefix
        << "[";

    for (
        std::size_t i = 0;
        i < P.coeff.size();
        ++i
    ) {
        if (i != 0) {
            std::cout << ",";
        }

        std::cout
            << P.coeff[i];
    }

    std::cout
        << "]\n";
}

int main() {
    std::cout
        << "START EXPERIMENT 304\n"
        << "COEFFICIENT-GCD FACTORIZATION\n"
        << "WHAT IS THE INTERNAL DIVISOR STRUCTURE OF THE PATH POLYNOMIALS?\n"
        << "\n";

    const std::vector<CaseData> cases = {
        {81077, 162749},
        {125017, 174259},
        {57107, 88261},
        {52757, 148457},
        {19483, 34123},
        {123001, 188291},
        {97987, 112583},
        {76129, 192113},
        {41257, 73643},
        {137239, 140419},
        {64187, 115319},
        {162527, 184087},
        {66179, 124753},
        {87943, 158047},
        {42683, 52529},
        {107473, 149711},
        {75853, 82759},
        {150131, 175267},
        {99859, 124769},
        {75337, 162229},
        {124471, 168043},
        {51563, 87683},
        {52237, 142123},
        {19069, 28751},
        {117043, 187637},
        {92203, 112067},
        {75617, 185869},
        {40823, 67843},
        {134369, 136709},
        {63667, 109211}
    };

    const std::vector<u64> bases = {
        2, 3, 5, 7, 11, 13
    };

    const std::vector<int> offsets = {
        -4, -3, -2, -1,
         0,
         1, 2, 3, 4
    };

    u64 total_polynomials = 0;

    u64 zero_polynomials = 0;
    u64 gcd_one = 0;
    u64 gcd_nontrivial = 0;

    u64 gcd_equals_p = 0;
    u64 gcd_equals_q = 0;
    u64 gcd_equals_N = 0;

    u64 proper_gcd = 0;

    u64 gcd_divides_p = 0;
    u64 gcd_divides_q = 0;
    u64 gcd_divides_N = 0;

    u64 p_only_root_count = 0;
    u64 q_only_root_count = 0;

    u64 first_nontrivial_printed = 0;

    /*
        Histogram for small prime factors of the coefficient gcd.
        This is deliberately limited to primes <= 1000 so we can
        quickly see repeated structural divisors.
    */
    std::vector<std::pair<u64, u64>> small_prime_frequency;

    auto increment_small_prime =
        [&small_prime_frequency](u64 prime) {
            for (auto &entry :
                 small_prime_frequency) {

                if (entry.first == prime) {
                    ++entry.second;
                    return;
                }
            }

            small_prime_frequency.push_back({
                prime,
                1
            });
        };

    for (
        std::size_t ci = 0;
        ci < cases.size();
        ++ci
    ) {
        const u64 p = cases[ci].p;
        const u64 q = cases[ci].q;
        const u64 N = p * q;

        u64 case_polynomials = 0;
        u64 case_zero = 0;
        u64 case_gcd_one = 0;
        u64 case_gcd_nontrivial = 0;
        u64 case_proper = 0;

        u64 case_p_only_roots = 0;
        u64 case_q_only_roots = 0;

        std::cout
            << "CASE "
            << ci
            << " p=" << p
            << " q=" << q
            << "\n";

        u64 root =
            integer_sqrt(N);

        for (u64 base : bases) {
            for (int offset : offsets) {

                i128 m_signed =
                    static_cast<i128>(root) +
                    static_cast<i128>(offset);

                if (m_signed <= 0) {
                    continue;
                }

                u64 m =
                    static_cast<u64>(
                        m_signed
                    );

                std::vector<u64> ns;

                ns.push_back(N);
                ns.push_back(m);
                ns.push_back(m + 1);

                if (
                    m != 0 &&
                    m <= N / m
                ) {
                    u64 D =
                        N - m * m;

                    if (D != 0) {
                        ns.push_back(D);
                    }
                }

                for (u64 n : ns) {
                    if (n == 0) {
                        continue;
                    }

                    i128 y =
                        static_cast<i128>(n) - 1;

                    Polynomial P =
                        build_polynomial(
                            m,
                            y,
                            base
                        );

                    ++case_polynomials;
                    ++total_polynomials;

                    u64 g =
                        coefficient_gcd(P);

                    if (g == 0) {
                        ++case_zero;
                        ++zero_polynomials;
                        continue;
                    }

                    if (g == 1) {
                        ++case_gcd_one;
                        ++gcd_one;
                    } else {
                        ++case_gcd_nontrivial;
                        ++gcd_nontrivial;
                    }

                    bool equals_p =
                        (g == p);

                    bool equals_q =
                        (g == q);

                    bool equals_N =
                        (g == N);

                    if (equals_p) {
                        ++gcd_equals_p;
                    }

                    if (equals_q) {
                        ++gcd_equals_q;
                    }

                    if (equals_N) {
                        ++gcd_equals_N;
                    }

                    bool divides_p =
                        (p % g == 0);

                    bool divides_q =
                        (q % g == 0);

                    bool divides_N =
                        (N % g == 0);

                    if (divides_p) {
                        ++gcd_divides_p;
                    }

                    if (divides_q) {
                        ++gcd_divides_q;
                    }

                    if (divides_N) {
                        ++gcd_divides_N;
                    }

                    bool is_proper =
                        (
                            g > 1 &&
                            g < N &&
                            g != p &&
                            g != q
                        );

                    if (is_proper) {
                        ++proper_gcd;
                        ++case_proper;
                    }

                    /*
                        Record small prime factors.
                    */
                    std::vector<Factor> fs =
                        factor_integer(g);

                    for (
                        const Factor &f : fs
                    ) {
                        if (f.prime <= 1000) {
                            increment_small_prime(
                                f.prime
                            );
                        }
                    }

                    /*
                        Only inspect roots for polynomials
                        that are not identically zero modulo
                        both primes.
                    */
                    bool div_p_coeff =
                        all_coeff_divisible(
                            P,
                            p
                        );

                    bool div_q_coeff =
                        all_coeff_divisible(
                            P,
                            q
                        );

                    if (
                        !div_p_coeff &&
                        !div_q_coeff
                    ) {
                        RootWitness Rp =
                            find_first_root(
                                P,
                                p,
                                g
                            );

                        RootWitness Rq =
                            find_first_root(
                                P,
                                q,
                                g
                            );

                        bool p_root =
                            Rp.found;

                        bool q_root =
                            Rq.found;

                        if (
                            p_root &&
                            !q_root
                        ) {
                            ++p_only_root_count;
                            ++case_p_only_roots;
                        }

                        if (
                            q_root &&
                            !p_root
                        ) {
                            ++q_only_root_count;
                            ++case_q_only_roots;
                        }

                        /*
                            Print all exceptional cases,
                            since there should be very few.
                        */
                        if (
                            (p_root && !q_root) ||
                            (q_root && !p_root)
                        ) {
                            std::cout
                                << "EXCEPTIONAL_ROOT\n"
                                << "type="
                                << (
                                    p_root
                                        ? "P_ONLY"
                                        : "Q_ONLY"
                                )
                                << "\n"
                                << "base="
                                << base
                                << " m="
                                << m
                                << " n="
                                << n
                                << "\n"
                                << "coefficient_gcd="
                                << g
                                << "\n"
                                << "factorization="
                                << factor_string(g)
                                << "\n"
                                << "root="
                                << (
                                    p_root
                                        ? Rp.root
                                        : Rq.root
                                )
                                << "\n"
                                << "gcd_with_N="
                                << gcd_u64(g, N)
                                << "\n";

                            print_coefficients(
                                P,
                                "coefficients="
                            );

                            /*
                                Stop dumping huge numbers of
                                exceptional cases if a bug
                                unexpectedly creates many.
                            */
                            if (
                                first_nontrivial_printed < 20
                            ) {
                                ++first_nontrivial_printed;
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "polynomials="
            << case_polynomials
            << "\n";

        std::cout
            << "zero_polynomials="
            << case_zero
            << "\n";

        std::cout
            << "gcd_one="
            << case_gcd_one
            << "\n";

        std::cout
            << "gcd_nontrivial="
            << case_gcd_nontrivial
            << "\n";

        std::cout
            << "proper_gcd="
            << case_proper
            << "\n";

        std::cout
            << "p_only_roots="
            << case_p_only_roots
            << "\n";

        std::cout
            << "q_only_roots="
            << case_q_only_roots
            << "\n\n";
    }

    std::sort(
        small_prime_frequency.begin(),
        small_prime_frequency.end(),
        [](
            const auto &a,
            const auto &b
        ) {
            if (a.second != b.second) {
                return a.second > b.second;
            }

            return a.first < b.first;
        }
    );

    std::cout
        << "============================\n"
        << "TOTAL\n";

    std::cout
        << "polynomials="
        << total_polynomials
        << "\n";

    std::cout
        << "zero_polynomials="
        << zero_polynomials
        << "\n";

    std::cout
        << "gcd_one="
        << gcd_one
        << "\n";

    std::cout
        << "gcd_nontrivial="
        << gcd_nontrivial
        << "\n";

    std::cout
        << "gcd_equals_p="
        << gcd_equals_p
        << "\n";

    std::cout
        << "gcd_equals_q="
        << gcd_equals_q
        << "\n";

    std::cout
        << "gcd_equals_N="
        << gcd_equals_N
        << "\n";

    std::cout
        << "proper_gcd="
        << proper_gcd
        << "\n";

    std::cout
        << "gcd_divides_p="
        << gcd_divides_p
        << "\n";

    std::cout
        << "gcd_divides_q="
        << gcd_divides_q
        << "\n";

    std::cout
        << "gcd_divides_N="
        << gcd_divides_N
        << "\n";

    std::cout
        << "p_only_root_count="
        << p_only_root_count
        << "\n";

    std::cout
        << "q_only_root_count="
        << q_only_root_count
        << "\n";

    std::cout
        << "\nSMALL PRIME FACTORS OF COEFFICIENT GCD\n";

    const std::size_t print_limit =
        std::min<std::size_t>(
            30,
            small_prime_frequency.size()
        );

    for (
        std::size_t i = 0;
        i < print_limit;
        ++i
    ) {
        std::cout
            << small_prime_frequency[i].first
            << " "
            << small_prime_frequency[i].second
            << "\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT 304\n";

    return 0;
}
