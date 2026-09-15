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

struct ReducedPolynomial {
    std::vector<u64> coeff;
};

struct NonzeroTerm {
    std::size_t degree;
    u64 coefficient;
};

struct Stats {
    u64 monomial = 0;
    u64 non_monomial = 0;

    u64 min_t_valuation = UINT64_MAX;
    u64 max_t_valuation = 0;

    u64 valuation_hist[32] = {};

    u64 support_size_sum = 0;
    u64 max_support_size = 0;

    u64 exceptional_printed = 0;
};

std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }

    bool negative = x < 0;

    u128 v =
        negative
            ? static_cast<u128>(-x)
            : static_cast<u128>(x);

    std::string s;

    while (v > 0) {
        int d =
            static_cast<int>(v % 10);

        s.push_back(
            static_cast<char>('0' + d)
        );

        v /= 10;
    }

    if (negative) {
        s.push_back('-');
    }

    std::reverse(
        s.begin(),
        s.end()
    );

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

void print_digits(
    u64 x,
    u64 base
) {
    std::vector<u64> d =
        base_digits(
            x,
            base
        );

    std::cout
        << "digits_lsf=[";

    for (
        std::size_t i = 0;
        i < d.size();
        ++i
    ) {
        if (i != 0) {
            std::cout << ",";
        }

        std::cout << d[i];
    }

    std::cout
        << "]\n";
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

    md.resize(
        L,
        0
    );

    yd.resize(
        L,
        0
    );

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
            answer +=
                suffix[pos + 1];
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
            << "ERROR: negative local derivative\n";

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
            << "ERROR: path did not reach target\n";

        std::exit(1);
    }

    return P;
}

u64 coefficient_gcd(
    const Polynomial &P
) {
    u64 g = 0;

    for (u64 c : P.coeff) {
        g = gcd_u64(g, c);
    }

    return g;
}

std::size_t t_valuation(
    const Polynomial &P
) {
    for (
        std::size_t i = 0;
        i < P.coeff.size();
        ++i
    ) {
        if (P.coeff[i] != 0) {
            return i;
        }
    }

    return P.coeff.size();
}

std::vector<NonzeroTerm> support(
    const Polynomial &P
) {
    std::vector<NonzeroTerm> result;

    for (
        std::size_t i = 0;
        i < P.coeff.size();
        ++i
    ) {
        if (P.coeff[i] != 0) {
            result.push_back({
                i,
                P.coeff[i]
            });
        }
    }

    return result;
}

ReducedPolynomial strip_all_t(
    const Polynomial &P,
    std::size_t k
) {
    ReducedPolynomial R;

    if (k >= P.coeff.size()) {
        return R;
    }

    R.coeff.resize(
        P.coeff.size() - k
    );

    for (
        std::size_t i = k;
        i < P.coeff.size();
        ++i
    ) {
        R.coeff[i - k] =
            P.coeff[i];
    }

    while (
        R.coeff.size() > 1 &&
        R.coeff.back() == 0
    ) {
        R.coeff.pop_back();
    }

    return R;
}

u64 reduced_support_size(
    const ReducedPolynomial &R
) {
    u64 count = 0;

    for (u64 c : R.coeff) {
        if (c != 0) {
            ++count;
        }
    }

    return count;
}

bool is_monomial(
    const Polynomial &P
) {
    u64 count = 0;

    for (u64 c : P.coeff) {
        if (c != 0) {
            ++count;
        }
    }

    return count == 1;
}

bool is_all_zero(
    const Polynomial &P
) {
    for (u64 c : P.coeff) {
        if (c != 0) {
            return false;
        }
    }

    return true;
}

void print_support(
    const Polynomial &P
) {
    std::cout
        << "support=[";

    bool first = true;

    for (
        std::size_t i = 0;
        i < P.coeff.size();
        ++i
    ) {
        if (P.coeff[i] == 0) {
            continue;
        }

        if (!first) {
            std::cout << ",";
        }

        first = false;

        std::cout
            << "("
            << i
            << ","
            << P.coeff[i]
            << ")";
    }

    std::cout
        << "]\n";
}

void accumulate(
    Stats &dst,
    const Stats &src
) {
    dst.monomial += src.monomial;
    dst.non_monomial += src.non_monomial;

    dst.min_t_valuation =
        std::min(
            dst.min_t_valuation,
            src.min_t_valuation
        );

    dst.max_t_valuation =
        std::max(
            dst.max_t_valuation,
            src.max_t_valuation
        );

    for (int i = 0; i < 32; ++i) {
        dst.valuation_hist[i] +=
            src.valuation_hist[i];
    }

    dst.support_size_sum +=
        src.support_size_sum;

    dst.max_support_size =
        std::max(
            dst.max_support_size,
            src.max_support_size
        );
}

int main() {
    std::cout
        << "START EXPERIMENT 307\n"
        << "FULL t-ADIC STRUCTURE OF PATH POLYNOMIALS\n"
        << "WHAT IS THE EXACT POWER OF t AND WHERE DO THE EXCEPTIONS LIVE?\n"
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

    Stats total;

    u64 total_polynomials = 0;
    u64 total_zero = 0;
    u64 total_nonzero = 0;

    u64 exceptional_count = 0;

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
        u64 case_nonzero = 0;
        u64 case_exceptional = 0;

        Stats case_stats;

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

                    if (is_all_zero(P)) {
                        ++case_zero;
                        ++total_zero;
                        continue;
                    }

                    ++case_nonzero;
                    ++total_nonzero;

                    std::vector<NonzeroTerm> terms =
                        support(P);

                    std::size_t k =
                        terms.front().degree;

                    u64 support_count =
                        static_cast<u64>(
                            terms.size()
                        );

                    bool monomial =
                        (support_count == 1);

                    Stats current;

                    if (monomial) {
                        current.monomial = 1;
                    } else {
                        current.non_monomial = 1;
                        ++case_exceptional;
                        ++exceptional_count;
                    }

                    current.min_t_valuation =
                        static_cast<u64>(k);

                    current.max_t_valuation =
                        static_cast<u64>(k);

                    if (k < 32) {
                        current.valuation_hist[k] = 1;
                    }

                    current.support_size_sum =
                        support_count;

                    current.max_support_size =
                        support_count;

                    accumulate(
                        case_stats,
                        current
                    );

                    accumulate(
                        total,
                        current
                    );

                    /*
                        Print every non-monomial polynomial.
                        There were only 18 in Experiment 306,
                        so this should remain small.
                    */
                    if (!monomial) {
                        std::cout
                            << "NON_MONOMIAL\n"
                            << "base="
                            << base
                            << " m="
                            << m
                            << " n="
                            << n
                            << "\n"
                            << "N="
                            << N
                            << "\n"
                            << "t_valuation="
                            << k
                            << "\n"
                            << "support_size="
                            << support_count
                            << "\n";

                        print_digits(
                            m,
                            base
                        );

                        print_support(P);

                        u64 g =
                            coefficient_gcd(P);

                        std::cout
                            << "content_gcd="
                            << g
                            << "\n";

                        ReducedPolynomial R =
                            strip_all_t(
                                P,
                                k
                            );

                        std::cout
                            << "after_stripping_t_degree="
                            << (
                                R.coeff.empty()
                                    ? 0
                                    : R.coeff.size() - 1
                            )
                            << "\n";

                        std::cout
                            << "Q=[";
                        
                        for (
                            std::size_t i = 0;
                            i < R.coeff.size();
                            ++i
                        ) {
                            if (i != 0) {
                                std::cout
                                    << ",";
                            }

                            std::cout
                                << R.coeff[i];
                        }

                        std::cout
                            << "]\n\n";
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
            << "nonzero_polynomials="
            << case_nonzero
            << "\n";

        std::cout
            << "monomial="
            << case_stats.monomial
            << "\n";

        std::cout
            << "non_monomial="
            << case_stats.non_monomial
            << "\n";

        std::cout
            << "min_t_valuation="
            << case_stats.min_t_valuation
            << "\n";

        std::cout
            << "max_t_valuation="
            << case_stats.max_t_valuation
            << "\n";

        std::cout
            << "max_support_size="
            << case_stats.max_support_size
            << "\n";

        std::cout
            << "\n";
    }

    std::cout
        << "============================\n"
        << "TOTAL\n";

    std::cout
        << "polynomials="
        << total_polynomials
        << "\n";

    std::cout
        << "zero_polynomials="
        << total_zero
        << "\n";

    std::cout
        << "nonzero_polynomials="
        << total_nonzero
        << "\n";

    std::cout
        << "monomial="
        << total.monomial
        << "\n";

    std::cout
        << "non_monomial="
        << total.non_monomial
        << "\n";

    std::cout
        << "exceptional_count="
        << exceptional_count
        << "\n";

    std::cout
        << "min_t_valuation="
        << total.min_t_valuation
        << "\n";

    std::cout
        << "max_t_valuation="
        << total.max_t_valuation
        << "\n";

    std::cout
        << "support_size_sum="
        << total.support_size_sum
        << "\n";

    std::cout
        << "max_support_size="
        << total.max_support_size
        << "\n";

    std::cout
        << "\nT-VALUATION HISTOGRAM\n";

    for (int k = 0; k < 32; ++k) {
        if (total.valuation_hist[k] != 0) {
            std::cout
                << k
                << " "
                << total.valuation_hist[k]
                << "\n";
        }
    }

    std::cout
        << "\nFINISHED EXPERIMENT 307\n";

    return 0;
}
