#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

using i128 = __int128_t;
using u64 = std::uint64_t;

struct Poly {
    std::vector<i128> c; // ascending powers
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
    return is_zero(a)
        ? -1
        : static_cast<int>(a.c.size()) - 1;
}

static Poly make_poly(
    std::initializer_list<i128> coeff
) {
    return trim(Poly{
        std::vector<i128>(coeff)
    });
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

static Poly poly_mod(
    Poly a,
    const Poly& b
) {
    a = trim(a);

    if (is_zero(b)) {
        return a;
    }

    const int db = degree(b);
    const i128 lead_b = b.c.back();

    while (!is_zero(a) && degree(a) >= db) {
        const int da = degree(a);

        /*
         * All cyclotomic divisors used here are monic.
         * Therefore exact integer division is valid.
         */
        const i128 factor =
            a.c.back() / lead_b;

        const int shift = da - db;

        for (int j = 0; j <= db; ++j) {
            a.c[j + shift] -=
                factor * b.c[j];
        }

        a = trim(a);
    }

    return a;
}

static bool exact_divisible(
    const Poly& f,
    const Poly& g
) {
    return is_zero(poly_mod(f, g));
}

static Poly multiply_x_mod(
    const Poly& a,
    const Poly& phi
) {
    const int d = degree(phi);

    std::vector<i128> out(
        static_cast<std::size_t>(d),
        0
    );

    for (int i = 0; i < static_cast<int>(a.c.size()); ++i) {
        const i128 v = a.c[i];

        if (v == 0) {
            continue;
        }

        if (i + 1 < d) {
            out[i + 1] += v;
        } else {
            /*
             * phi is monic:
             *
             * x^d = -(phi_0 + ... + phi_{d-1}x^{d-1})
             */
            for (int j = 0; j < d; ++j) {
                out[j] -=
                    v * phi.c[j];
            }
        }
    }

    return trim(Poly{std::move(out)});
}

static i128 determinant_bareiss(
    std::vector<std::vector<i128>> a
) {
    const int n =
        static_cast<int>(a.size());

    if (n == 0) {
        return 1;
    }

    if (n == 1) {
        return a[0][0];
    }

    i128 prev_pivot = 1;
    int sign = 1;

    for (int k = 0; k < n - 1; ++k) {
        int pivot = k;

        while (
            pivot < n &&
            a[pivot][k] == 0
        ) {
            ++pivot;
        }

        if (pivot == n) {
            return 0;
        }

        if (pivot != k) {
            std::swap(a[pivot], a[k]);
            sign = -sign;
        }

        const i128 pivot_value =
            a[k][k];

        for (int i = k + 1; i < n; ++i) {
            for (int j = k + 1; j < n; ++j) {
                const i128 numerator =
                    a[i][j] * pivot_value
                    - a[i][k] * a[k][j];

                a[i][j] =
                    numerator / prev_pivot;
            }
        }

        for (int i = k + 1; i < n; ++i) {
            a[i][k] = 0;
        }

        for (int j = k + 1; j < n; ++j) {
            a[k][j] = a[k][j];
        }

        prev_pivot = pivot_value;
    }

    const i128 det =
        a[n - 1][n - 1];

    return sign < 0 ? -det : det;
}

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
        std::vector<i128>(
            size,
            0
        )
    );

    /*
     * First n rows: shifted copies of A.
     */
    for (int row = 0; row < n; ++row) {
        for (int j = 0; j <= m; ++j) {
            S[row][row + j] =
                A.c[j];
        }
    }

    /*
     * Next m rows: shifted copies of B.
     */
    for (int row = 0; row < m; ++row) {
        for (int j = 0; j <= n; ++j) {
            S[n + row][row + j] =
                B.c[j];
        }
    }

    return determinant_bareiss(S);
}

static i128 resultant_multiplication_map(
    const Poly& f,
    const Poly& phi
) {
    const int d = degree(phi);

    if (d < 0) {
        return 0;
    }

    /*
     * r = f mod phi.
     */
    const Poly r =
        poly_mod(f, phi);

    std::vector<std::vector<i128>> M(
        d,
        std::vector<i128>(
            d,
            0
        )
    );

    /*
     * Columns are multiplication by
     *
     *   1, x, x^2, ..., x^(d-1)
     *
     * inside Z[x]/(phi).
     */
    Poly basis =
        make_poly({1});

    for (int col = 0; col < d; ++col) {
        const Poly product =
            poly_mod(
                poly_mul(r, basis),
                phi
            );

        for (int row = 0; row < d; ++row) {
            if (row <
                static_cast<int>(
                    product.c.size()
                )) {
                M[row][col] =
                    product.c[row];
            }
        }

        basis =
            multiply_x_mod(
                basis,
                phi
            );
    }

    return determinant_bareiss(M);
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

    return trim(Poly{
        std::move(c)
    });
}

static bool same_up_to_sign(
    i128 a,
    i128 b
) {
    return a == b || a == -b;
}

static std::string to_string_i128(
    i128 x
) {
    if (x == 0) {
        return "0";
    }

    bool negative = false;

    if (x < 0) {
        negative = true;
        x = -x;
    }

    std::string s;

    while (x > 0) {
        s.push_back(
            static_cast<char>(
                '0' +
                static_cast<int>(
                    x % 10
                )
            )
        );

        x /= 10;
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

int main() {
    std::cout << "START EXPERIMENT 320\n";
    std::cout << "FAST RESULTANT IMPLEMENTATION AUDIT\n";
    std::cout << "MULTIPLICATION-MAP VS BAREISS SYLVESTER RESULTANTS\n\n";

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
     * --------------------------------------------------
     * TEST 1: x^n - 1
     * --------------------------------------------------
     */
    u64 known_tests = 0;
    u64 known_failures = 0;

    for (const auto& phi : phis) {
        for (int n = 1; n <= 20; ++n) {
            const Poly f =
                x_pow_minus_one(n);

            const i128 r_map =
                resultant_multiplication_map(
                    f,
                    phi.poly
                );

            const i128 r_syl =
                resultant_sylvester(
                    f,
                    phi.poly
                );

            ++known_tests;

            if (!same_up_to_sign(
                    r_map,
                    r_syl)) {

                ++known_failures;

                std::cout
                    << "KNOWN_RESULTANT_FAIL"
                    << " phi=" << phi.k
                    << " n=" << n
                    << " map="
                    << to_string_i128(r_map)
                    << " syl="
                    << to_string_i128(r_syl)
                    << "\n";
            }

            /*
             * Phi_k | x^n - 1 iff k | n.
             */
            const bool expected_zero =
                (n % phi.k == 0);

            const bool actual_zero =
                (r_syl == 0);

            if (expected_zero != actual_zero) {
                ++known_failures;

                std::cout
                    << "KNOWN_ZERO_FAIL"
                    << " phi=" << phi.k
                    << " n=" << n
                    << " expected="
                    << expected_zero
                    << " actual="
                    << actual_zero
                    << "\n";
            }
        }
    }

    std::cout
        << "KNOWN TESTS\n"
        << "tests=" << known_tests << "\n"
        << "failures=" << known_failures << "\n\n";

    /*
     * --------------------------------------------------
     * TEST 2: deliberately constructed divisible cases
     * --------------------------------------------------
     */
    u64 constructed_tests = 0;
    u64 constructed_failures = 0;

    for (const auto& phi : phis) {
        for (int seed = 1; seed <= 30; ++seed) {
            std::vector<i128> qcoeff(
                static_cast<std::size_t>(
                    2 + (seed % 5)
                ),
                0
            );

            for (std::size_t i = 0;
                 i < qcoeff.size();
                 ++i) {

                qcoeff[i] =
                    static_cast<i128>(
                        static_cast<int>(
                            (seed + 7 * i) % 11
                        ) - 5
                    );
            }

            if (qcoeff.back() == 0) {
                qcoeff.back() = 1;
            }

            const Poly quotient =
                trim(Poly{
                    std::move(qcoeff)
                });

            const Poly f =
                poly_mul(
                    phi.poly,
                    quotient
                );

            const i128 r_map =
                resultant_multiplication_map(
                    f,
                    phi.poly
                );

            const i128 r_syl =
                resultant_sylvester(
                    f,
                    phi.poly
                );

            ++constructed_tests;

            if (r_map != 0 ||
                r_syl != 0 ||
                !exact_divisible(
                    f,
                    phi.poly
                )) {

                ++constructed_failures;

                std::cout
                    << "CONSTRUCTED_FAIL"
                    << " phi=" << phi.k
                    << " seed=" << seed
                    << " map="
                    << to_string_i128(r_map)
                    << " syl="
                    << to_string_i128(r_syl)
                    << "\n";
            }
        }
    }

    std::cout
        << "CONSTRUCTED DIVISIBLE TESTS\n"
        << "tests=" << constructed_tests << "\n"
        << "failures=" << constructed_failures
        << "\n\n";

    /*
     * --------------------------------------------------
     * TEST 3: small non-divisible cases
     * --------------------------------------------------
     */
    u64 nondiv_tests = 0;
    u64 nondiv_failures = 0;

    for (const auto& phi : phis) {
        for (int seed = 1; seed <= 30; ++seed) {
            const int d =
                1 + (seed % 7);

            std::vector<i128> coeff(
                static_cast<std::size_t>(
                    d + 1
                ),
                0
            );

            for (int i = 0; i <= d; ++i) {
                coeff[i] =
                    static_cast<i128>(
                        static_cast<int>(
                            (seed * 3 + i * 5)
                            % 13
                        ) - 6
                    );
            }

            if (coeff.back() == 0) {
                coeff.back() = 1;
            }

            const Poly f =
                trim(Poly{
                    std::move(coeff)
                });

            const bool divisible =
                exact_divisible(
                    f,
                    phi.poly
                );

            const i128 r_map =
                resultant_multiplication_map(
                    f,
                    phi.poly
                );

            const i128 r_syl =
                resultant_sylvester(
                    f,
                    phi.poly
                );

            ++nondiv_tests;

            if (!same_up_to_sign(
                    r_map,
                    r_syl) ||
                divisible != (r_syl == 0)) {

                ++nondiv_failures;

                std::cout
                    << "NONDIV_FAIL"
                    << " phi=" << phi.k
                    << " seed=" << seed
                    << " divisible="
                    << divisible
                    << " map="
                    << to_string_i128(
                           r_map)
                    << " syl="
                    << to_string_i128(
                           r_syl)
                    << "\n";
            }
        }
    }

    std::cout
        << "SMALL GENERAL TESTS\n"
        << "tests=" << nondiv_tests << "\n"
        << "failures=" << nondiv_failures
        << "\n\n";

    /*
     * --------------------------------------------------
     * TOTAL
     * --------------------------------------------------
     */
    const u64 total_failures =
        known_failures +
        constructed_failures +
        nondiv_failures;

    std::cout
        << "============================\n"
        << "TOTAL\n"
        << "known_tests="
        << known_tests
        << "\n"
        << "constructed_tests="
        << constructed_tests
        << "\n"
        << "general_tests="
        << nondiv_tests
        << "\n"
        << "total_failures="
        << total_failures
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 320\n";

    return 0;
}
