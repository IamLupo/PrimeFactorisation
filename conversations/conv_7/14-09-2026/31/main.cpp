#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using i128 = __int128_t;
using u64 = std::uint64_t;

struct Poly {
    std::vector<i128> c; // ascending
};

struct Cyclotomic {
    int k;
    Poly poly;
};

static Poly trim(Poly a) {
    while (a.c.size() > 1 && a.c.back() == 0) {
        a.c.pop_back();
    }
    if (a.c.empty()) {
        a.c.push_back(0);
    }
    return a;
}

static bool is_zero(const Poly& a) {
    return a.c.size() == 1 && a.c[0] == 0;
}

static int degree(const Poly& a) {
    if (is_zero(a)) {
        return -1;
    }
    return static_cast<int>(a.c.size()) - 1;
}

static Poly make_poly(std::initializer_list<i128> x) {
    return trim(Poly{std::vector<i128>(x)});
}

static Cyclotomic make_cyclotomic(int k) {
    switch (k) {
        case 2:
            return {2, make_poly({1, 1})};

        case 3:
            return {3, make_poly({1, 1, 1})};

        case 4:
            return {4, make_poly({1, 0, 1})};

        case 5:
            return {5, make_poly({1, 1, 1, 1, 1})};

        case 6:
            return {6, make_poly({1, -1, 1})};

        case 8:
            return {8, make_poly({1, 0, 0, 0, 1})};

        case 10:
            return {10, make_poly({1, -1, 1, -1, 1})};

        case 12:
            return {12, make_poly({1, 0, -1, 0, 1})};

        default:
            return {0, make_poly({0})};
    }
}

static Poly poly_mod(
    Poly a,
    const Poly& b
) {
    a = trim(a);

    if (is_zero(b)) {
        return a;
    }

    const int db = degree(b);

    while (!is_zero(a) && degree(a) >= db) {
        const int da = degree(a);
        const i128 lead = a.c.back() / b.c.back();
        const int shift = da - db;

        for (int j = 0; j <= db; ++j) {
            a.c[j + shift] -= lead * b.c[j];
        }

        a = trim(a);
    }

    return a;
}

static bool exact_divisible(
    const Poly& f,
    const Poly& g
) {
    if (is_zero(g)) {
        return false;
    }

    Poly r = poly_mod(f, g);
    return is_zero(r);
}

static Poly poly_mul(
    const Poly& a,
    const Poly& b
) {
    std::vector<i128> c(
        a.c.size() + b.c.size() - 1,
        0
    );

    for (std::size_t i = 0; i < a.c.size(); ++i) {
        for (std::size_t j = 0; j < b.c.size(); ++j) {
            c[i + j] += a.c[i] * b.c[j];
        }
    }

    return trim(Poly{std::move(c)});
}

static Poly poly_shift(
    const Poly& a,
    int shift
) {
    if (is_zero(a)) {
        return a;
    }

    std::vector<i128> c(
        a.c.size() + shift,
        0
    );

    for (std::size_t i = 0; i < a.c.size(); ++i) {
        c[i + shift] = a.c[i];
    }

    return Poly{std::move(c)};
}

static Poly poly_sub(
    const Poly& a,
    const Poly& b
) {
    const std::size_t n =
        std::max(a.c.size(), b.c.size());

    std::vector<i128> c(n, 0);

    for (std::size_t i = 0; i < n; ++i) {
        const i128 x =
            i < a.c.size() ? a.c[i] : 0;

        const i128 y =
            i < b.c.size() ? b.c[i] : 0;

        c[i] = x - y;
    }

    return trim(Poly{std::move(c)});
}

/*
 * Exact determinant by permutation.
 * Degrees here are tiny, so this is only used for
 * the independent multiplication-map calculation.
 */
static i128 determinant(
    const std::vector<std::vector<i128>>& A
) {
    const int n = static_cast<int>(A.size());

    if (n == 0) {
        return 1;
    }

    std::vector<int> p(n);

    for (int i = 0; i < n; ++i) {
        p[i] = i;
    }

    i128 det = 0;

    do {
        i128 term = 1;
        int inversions = 0;

        for (int i = 0; i < n; ++i) {
            term *= A[i][p[i]];

            for (int j = i + 1; j < n; ++j) {
                if (p[i] > p[j]) {
                    ++inversions;
                }
            }
        }

        if ((inversions & 1) == 0) {
            det += term;
        } else {
            det -= term;
        }

    } while (std::next_permutation(p.begin(), p.end()));

    return det;
}

/*
 * x * a(x) modulo monic phi.
 */
static Poly multiply_x_mod(
    const Poly& a,
    const Poly& phi
) {
    const int d = degree(phi);

    std::vector<i128> out(
        static_cast<std::size_t>(d),
        0
    );

    for (int i = 0; i < degree(a) + 1; ++i) {
        const i128 v = a.c[i];

        if (v == 0) {
            continue;
        }

        if (i + 1 < d) {
            out[i + 1] += v;
        } else {
            /*
             * Since phi is monic:
             *
             * x^d = -(phi_0 + ... + phi_{d-1}x^{d-1})
             */
            for (int j = 0; j < d; ++j) {
                out[j] -= v * phi.c[j];
            }
        }
    }

    return trim(Poly{std::move(out)});
}

static Poly reduce_mod_phi(
    const Poly& f,
    const Poly& phi
) {
    return poly_mod(f, phi);
}

/*
 * Version equivalent to the multiplication-map approach
 * used in Experiment 318, but implemented from Poly objects.
 */
static i128 resultant_multiplication_map(
    const Poly& f,
    const Poly& phi
) {
    const int d = degree(phi);

    if (d < 0) {
        return 0;
    }

    const Poly r =
        reduce_mod_phi(f, phi);

    std::vector<std::vector<i128>> M(
        d,
        std::vector<i128>(d, 0)
    );

    Poly basis = make_poly({1});

    for (int col = 0; col < d; ++col) {
        Poly product =
            poly_mod(
                poly_mul(r, basis),
                phi
            );

        for (int row = 0; row < d; ++row) {
            if (row < static_cast<int>(product.c.size())) {
                M[row][col] = product.c[row];
            }
        }

        basis =
            multiply_x_mod(basis, phi);
    }

    return determinant(M);
}

/*
 * Direct resultant via the Sylvester matrix.
 *
 * For modest degrees this is an excellent independent
 * implementation and is deliberately unrelated to the
 * multiplication-map construction.
 */
static i128 resultant_sylvester(
    const Poly& A,
    const Poly& B
) {
    const int m = degree(A);
    const int n = degree(B);

    if (m < 0 || n < 0) {
        return 0;
    }

    if (m == 0) {
        /*
         * Res(c,B)=c^deg(B)
         */
        i128 result = 1;

        for (int i = 0; i < n; ++i) {
            result *= A.c[0];
        }

        return result;
    }

    if (n == 0) {
        i128 result = 1;

        for (int i = 0; i < m; ++i) {
            result *= B.c[0];
        }

        return result;
    }

    const int size = m + n;

    std::vector<std::vector<i128>> S(
        size,
        std::vector<i128>(size, 0)
    );

    /*
     * First n rows: A shifted.
     */
    for (int row = 0; row < n; ++row) {
        for (int j = 0; j <= m; ++j) {
            S[row][row + j] = A.c[j];
        }
    }

    /*
     * Next m rows: B shifted.
     */
    for (int row = 0; row < m; ++row) {
        for (int j = 0; j <= n; ++j) {
            S[n + row][row + j] = B.c[j];
        }
    }

    return determinant(S);
}

static Poly x_pow_minus_one(
    int n
) {
    std::vector<i128> c(
        static_cast<std::size_t>(n + 1),
        0
    );

    c[0] = -1;
    c[n] = 1;

    return trim(Poly{std::move(c)});
}

static std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }

    bool neg = x < 0;

    if (neg) {
        x = -x;
    }

    std::string s;

    while (x > 0) {
        s.push_back(
            static_cast<char>(
                '0' + static_cast<int>(x % 10)
            )
        );

        x /= 10;
    }

    if (neg) {
        s.push_back('-');
    }

    std::reverse(s.begin(), s.end());

    return s;
}

static bool same_up_to_sign(
    i128 a,
    i128 b
) {
    return a == b || a == -b;
}

int main() {
    std::cout << "START EXPERIMENT 320\n";
    std::cout << "RESULTANT IMPLEMENTATION AUDIT\n";
    std::cout << "DOES THE MULTIPLICATION-MAP RESULTANT AGREE WITH AN INDEPENDENT SYLVESTER RESULTANT?\n\n";

    const std::array<Cyclotomic, 8> phis = {
        make_cyclotomic(2),
        make_cyclotomic(3),
        make_cyclotomic(4),
        make_cyclotomic(5),
        make_cyclotomic(6),
        make_cyclotomic(8),
        make_cyclotomic(10),
        make_cyclotomic(12)
    };

    /*
     * 1. Basic known-polynomial tests.
     */
    u64 known_tests = 0;
    u64 known_failures = 0;

    std::cout << "KNOWN RESULTANT TESTS\n";

    for (const auto& phi : phis) {
        /*
         * If Phi_k | x^n-1 exactly when k | n,
         * then Res(x^n-1,Phi_k)=0.
         */
        for (int n = 1; n <= 24; ++n) {
            const Poly f =
                x_pow_minus_one(n);

            const i128 r1 =
                resultant_multiplication_map(
                    f,
                    phi.poly
                );

            const i128 r2 =
                resultant_sylvester(
                    f,
                    phi.poly
                );

            ++known_tests;

            /*
             * The two implementations must agree up
             * to the conventional sign/order.
             */
            if (!same_up_to_sign(r1, r2)) {
                ++known_failures;

                std::cout
                    << "KNOWN_FAIL"
                    << " phi=" << phi.k
                    << " n=" << n
                    << " map=" << to_string_i128(r1)
                    << " sylvester="
                    << to_string_i128(r2)
                    << "\n";
            }

            const bool expected_zero =
                (n % phi.k == 0);

            const bool actual_zero =
                (r2 == 0);

            if (expected_zero != actual_zero) {
                ++known_failures;

                std::cout
                    << "KNOWN_ZERO_FAIL"
                    << " phi=" << phi.k
                    << " n=" << n
                    << " expected_zero="
                    << expected_zero
                    << " actual_zero="
                    << actual_zero
                    << "\n";
            }
        }
    }

    std::cout
        << "known_tests="
        << known_tests
        << "\n"
        << "known_failures="
        << known_failures
        << "\n\n";

    /*
     * 2. Random small polynomial agreement.
     */
    u64 random_tests = 0;
    u64 random_failures = 0;

    std::uint64_t state =
        0x9E3779B97F4A7C15ULL;

    auto rng = [&]() {
        state ^= state << 7;
        state ^= state >> 9;
        state ^= state << 8;
        return state;
    };

    for (int test = 0; test < 2000; ++test) {
        const auto& phi =
            phis[rng() % phis.size()];

        const int degree_f =
            1 + static_cast<int>(rng() % 10);

        std::vector<i128> c(
            static_cast<std::size_t>(degree_f + 1),
            0
        );

        for (int i = 0; i <= degree_f; ++i) {
            c[i] =
                static_cast<i128>(
                    static_cast<int>(rng() % 11) - 5
                );
        }

        if (c.back() == 0) {
            c.back() = 1;
        }

        const Poly f = trim(Poly{std::move(c)});

        const i128 r1 =
            resultant_multiplication_map(
                f,
                phi.poly
            );

        const i128 r2 =
            resultant_sylvester(
                f,
                phi.poly
            );

        ++random_tests;

        if (!same_up_to_sign(r1, r2)) {
            ++random_failures;

            if (random_failures <= 20) {
                std::cout
                    << "RANDOM_FAIL"
                    << " phi=" << phi.k
                    << " map="
                    << to_string_i128(r1)
                    << " sylvester="
                    << to_string_i128(r2)
                    << "\n";
            }
        }
    }

    std::cout
        << "random_tests="
        << random_tests
        << "\n"
        << "random_failures="
        << random_failures
        << "\n\n";

    /*
     * 3. Direct divisibility/resultant consistency.
     */
    u64 consistency_tests = 0;
    u64 consistency_failures = 0;

    for (const auto& phi : phis) {
        for (int n = 1; n <= 30; ++n) {
            Poly f = x_pow_minus_one(n);

            const bool divisible =
                exact_divisible(
                    f,
                    phi.poly
                );

            const i128 r =
                resultant_sylvester(
                    f,
                    phi.poly
                );

            const bool zero =
                (r == 0);

            ++consistency_tests;

            if (divisible != zero) {
                ++consistency_failures;

                std::cout
                    << "CONSISTENCY_FAIL"
                    << " phi=" << phi.k
                    << " n=" << n
                    << " divisible="
                    << divisible
                    << " resultant_zero="
                    << zero
                    << "\n";
            }
        }
    }

    std::cout
        << "consistency_tests="
        << consistency_tests
        << "\n"
        << "consistency_failures="
        << consistency_failures
        << "\n\n";

    /*
     * 4. Small structural samples corresponding to 319.
     *
     * These are deliberately independent from the
     * semiprime set: we just generate arbitrary polynomials
     * and verify that zero resultant <=> exact divisibility
     * for these cyclotomics.
     */
    u64 structural_tests = 0;
    u64 structural_failures = 0;

    for (const auto& phi : phis) {
        for (int seed = 0; seed < 100; ++seed) {
            const int d =
                1 + (seed % 8);

            std::vector<i128> c(
                static_cast<std::size_t>(d + 1),
                0
            );

            for (int i = 0; i <= d; ++i) {
                c[i] =
                    static_cast<i128>(
                        ((seed + 3 * i) % 9) - 4
                    );
            }

            if (c.back() == 0) {
                c.back() = 1;
            }

            Poly f = trim(Poly{std::move(c)});

            /*
             * Force every 10th sample to contain Phi
             * exactly, giving us a guaranteed zero case.
             */
            if (seed % 10 == 0) {
                f = poly_mul(
                    phi.poly,
                    f
                );
            }

            const bool divisible =
                exact_divisible(
                    f,
                    phi.poly
                );

            const i128 r1 =
                resultant_multiplication_map(
                    f,
                    phi.poly
                );

            const i128 r2 =
                resultant_sylvester(
                    f,
                    phi.poly
                );

            ++structural_tests;

            if (!same_up_to_sign(r1, r2) ||
                divisible != (r2 == 0)) {

                ++structural_failures;

                if (structural_failures <= 20) {
                    std::cout
                        << "STRUCTURAL_FAIL"
                        << " phi=" << phi.k
                        << " seed=" << seed
                        << " divisible="
                        << divisible
                        << " map="
                        << to_string_i128(r1)
                        << " sylvester="
                        << to_string_i128(r2)
                        << "\n";
                }
            }
        }
    }

    std::cout
        << "structural_tests="
        << structural_tests
        << "\n"
        << "structural_failures="
        << structural_failures
        << "\n";

    std::cout
        << "\n============================\n"
        << "TOTAL\n"
        << "known_failures="
        << known_failures
        << "\n"
        << "random_failures="
        << random_failures
        << "\n"
        << "consistency_failures="
        << consistency_failures
        << "\n"
        << "structural_failures="
        << structural_failures
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 320\n";

    return 0;
}
