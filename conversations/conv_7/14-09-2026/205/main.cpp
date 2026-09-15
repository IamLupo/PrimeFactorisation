#include <algorithm>
#include <chrono>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

using namespace std;

static const int PRIME_LIMIT = 100000;
static const int J_LIMIT = 97;

struct Factorization
{
    vector<int> primes;
    vector<int> multiplicities;
};

struct CanonicalSequence
{
    int index = 0;
    vector<int> values;
};

vector<int> build_primes(int limit);
vector<int> build_spf(int limit);

bool is_prime(int n);
vector<int> factorize_distinct(int n);

Factorization factorize(
    int n,
    const vector<int>& spf
);

int mod_pow(
    int a,
    int e,
    int mod
);

bool is_primitive_root(
    int g,
    int p
);

int primitive_root_prime(int p);

vector<int> build_discrete_log_table(
    int g,
    int p
);

vector<int> all_divisors(
    const Factorization& f
);

bool has_proper_one_residue_divisor(
    int n,
    int j,
    const vector<int>& divisors
);

bool has_proper_zero_sum_subset(
    const vector<int>& exponents,
    int modulus
);

CanonicalSequence canonicalize(
    const vector<int>& exponents,
    int modulus
);

bool has_subset_sum(
    const vector<int>& values,
    int target
);

int distinct_count(
    const vector<int>& values
);

void print_sequence(
    const vector<int>& values
);

void print_factorization(
    const Factorization& f
);

void main_experiment();

int main()
{
    cout << "START EXPERIMENT 500\n";
    main_experiment();
    cout << "FINISHED EXPERIMENT 500\n";
    return 0;
}

vector<int> build_primes(int limit)
{
    vector<bool> composite(limit + 1, false);
    vector<int> primes;

    for (int i = 2; i <= limit; ++i)
    {
        if (composite[i])
        {
            continue;
        }

        primes.push_back(i);

        if (1LL * i * i <= limit)
        {
            for (int j = i * i; j <= limit; j += i)
            {
                composite[j] = true;
            }
        }
    }

    return primes;
}

vector<int> build_spf(int limit)
{
    vector<int> spf(limit + 1, 0);

    for (int i = 2; i <= limit; ++i)
    {
        if (spf[i] != 0)
        {
            continue;
        }

        spf[i] = i;

        if (1LL * i * i <= limit)
        {
            for (int j = i * i; j <= limit; j += i)
            {
                if (spf[j] == 0)
                {
                    spf[j] = i;
                }
            }
        }
    }

    return spf;
}

bool is_prime(int n)
{
    if (n < 2)
    {
        return false;
    }

    if (n % 2 == 0)
    {
        return n == 2;
    }

    for (int d = 3; 1LL * d * d <= n; d += 2)
    {
        if (n % d == 0)
        {
            return false;
        }
    }

    return true;
}

vector<int> factorize_distinct(int n)
{
    vector<int> factors;

    for (int p = 2; 1LL * p * p <= n; ++p)
    {
        if (n % p != 0)
        {
            continue;
        }

        factors.push_back(p);

        while (n % p == 0)
        {
            n /= p;
        }
    }

    if (n > 1)
    {
        factors.push_back(n);
    }

    return factors;
}

Factorization factorize(
    int n,
    const vector<int>& spf
)
{
    Factorization result;

    while (n > 1)
    {
        const int p = spf[n];

        int exponent = 0;

        while (n % p == 0)
        {
            n /= p;
            ++exponent;
        }

        result.primes.push_back(p);
        result.multiplicities.push_back(exponent);
    }

    return result;
}

int mod_pow(
    int a,
    int e,
    int mod
)
{
    long long result = 1;
    long long base = a % mod;

    while (e > 0)
    {
        if (e & 1)
        {
            result = (result * base) % mod;
        }

        base = (base * base) % mod;
        e >>= 1;
    }

    return static_cast<int>(result);
}

bool is_primitive_root(
    int g,
    int p
)
{
    if (p == 2)
    {
        return g == 1;
    }

    const int phi = p - 1;

    for (int q : factorize_distinct(phi))
    {
        if (mod_pow(g, phi / q, p) == 1)
        {
            return false;
        }
    }

    return true;
}

int primitive_root_prime(int p)
{
    if (p == 2)
    {
        return 1;
    }

    for (int g = 2; g < p; ++g)
    {
        if (is_primitive_root(g, p))
        {
            return g;
        }
    }

    return -1;
}

vector<int> build_discrete_log_table(
    int g,
    int p
)
{
    vector<int> logs(p, -1);

    if (p == 2)
    {
        logs[1] = 0;
        return logs;
    }

    long long value = 1;

    for (int e = 0; e < p - 1; ++e)
    {
        logs[static_cast<int>(value)] = e;
        value = (value * g) % p;
    }

    return logs;
}

vector<int> all_divisors(
    const Factorization& f
)
{
    vector<int> divisors = {1};

    for (size_t i = 0; i < f.primes.size(); ++i)
    {
        const int p = f.primes[i];
        const int exponent = f.multiplicities[i];

        vector<int> next;

        int power = 1;

        for (int e = 0; e <= exponent; ++e)
        {
            for (int d : divisors)
            {
                next.push_back(d * power);
            }

            power *= p;
        }

        divisors.swap(next);
    }

    sort(divisors.begin(), divisors.end());

    return divisors;
}

bool has_proper_one_residue_divisor(
    int n,
    int j,
    const vector<int>& divisors
)
{
    for (int d : divisors)
    {
        if (d <= 1 || d >= n)
        {
            continue;
        }

        if (d % j == 1)
        {
            return true;
        }
    }

    return false;
}

bool has_proper_zero_sum_subset(
    const vector<int>& exponents,
    int modulus
)
{
    const int omega =
        static_cast<int>(exponents.size());

    if (omega <= 1)
    {
        return false;
    }

    vector<vector<bool>> dp(
        omega + 1,
        vector<bool>(modulus, false)
    );

    dp[0][0] = true;

    for (int x : exponents)
    {
        vector<vector<bool>> next = dp;

        for (int count = 0; count < omega; ++count)
        {
            for (int sum = 0; sum < modulus; ++sum)
            {
                if (!dp[count][sum])
                {
                    continue;
                }

                next[count + 1]
                    [(sum + x) % modulus] = true;
            }
        }

        dp.swap(next);
    }

    for (int count = 1; count < omega; ++count)
    {
        if (dp[count][0])
        {
            return true;
        }
    }

    return false;
}

bool has_subset_sum(
    const vector<int>& values,
    int target
)
{
    vector<bool> dp(target + 1, false);
    dp[0] = true;

    for (int x : values)
    {
        if (x > target)
        {
            continue;
        }

        for (int s = target; s >= x; --s)
        {
            if (dp[s - x])
            {
                dp[s] = true;
            }
        }
    }

    return dp[target];
}

CanonicalSequence canonicalize(
    const vector<int>& exponents,
    int modulus
)
{
    CanonicalSequence best;

    best.index = 1000000000;

    for (int u = 1; u < modulus; ++u)
    {
        if (gcd(u, modulus) != 1)
        {
            continue;
        }

        vector<int> values;

        values.reserve(exponents.size());

        int sum = 0;

        for (int e : exponents)
        {
            const int x =
                (u * e) % modulus;

            if (x == 0)
            {
                values.clear();
                break;
            }

            values.push_back(x);
            sum += x;
        }

        if (values.empty())
        {
            continue;
        }

        if (sum % modulus != 0)
        {
            continue;
        }

        const int index =
            sum / modulus;

        sort(values.begin(), values.end());

        /*
         * Primary criterion: minimum index.
         *
         * Secondary criterion: lexicographically smallest
         * sorted representative, giving us a reproducible
         * canonical form.
         */
        if (index < best.index ||
            (index == best.index &&
             values < best.values))
        {
            best.index = index;
            best.values = values;
        }
    }

    return best;
}

int distinct_count(
    const vector<int>& values
)
{
    vector<int> copy = values;

    sort(copy.begin(), copy.end());

    copy.erase(
        unique(copy.begin(), copy.end()),
        copy.end()
    );

    return static_cast<int>(copy.size());
}

void print_sequence(
    const vector<int>& values
)
{
    cout << "[";

    for (size_t i = 0; i < values.size(); ++i)
    {
        if (i != 0)
        {
            cout << ",";
        }

        cout << values[i];
    }

    cout << "]";
}

void print_factorization(
    const Factorization& f
)
{
    for (size_t i = 0; i < f.primes.size(); ++i)
    {
        if (i != 0)
        {
            cout << "*";
        }

        cout << f.primes[i];

        if (f.multiplicities[i] > 1)
        {
            cout << "^" << f.multiplicities[i];
        }
    }
}

void main_experiment()
{
    const auto start =
        chrono::high_resolution_clock::now();

    const vector<int> primes =
        build_primes(PRIME_LIMIT);

    const int max_n =
        J_LIMIT * PRIME_LIMIT + 1;

    const vector<int> spf =
        build_spf(max_n);

    cout << "PRIME_COUNT="
         << primes.size()
         << "\n";

    cout << "J_LIMIT="
         << J_LIMIT
         << "\n";

    cout << "MAX_JP_PLUS_1="
         << max_n
         << "\n";

    cout << "SPF_READY=1\n";

    long long total_failures = 0;

    long long index_one = 0;
    long long index_two = 0;
    long long index_three = 0;
    long long index_four_plus = 0;

    long long index_two_subset_m =
        0;

    long long index_two_subset_m_omega =
        0;

    long long index_two_no_subset_m =
        0;

    long long total_index = 0;
    long long total_omega = 0;

    for (int j = 2; j <= J_LIMIT; ++j)
    {
        if (!is_prime(j))
        {
            continue;
        }

        const int modulus = j - 1;

        const int primitive_root =
            primitive_root_prime(j);

        if (primitive_root < 0)
        {
            cout << "PRIMITIVE_ROOT_FAILURE J="
                 << j
                 << "\n";

            return;
        }

        const vector<int> log_table =
            build_discrete_log_table(
                primitive_root,
                j
            );

        long long j_failures = 0;
        long long j_index_two = 0;

        for (int p : primes)
        {
            const long long n64 =
                1LL * j * p + 1;

            if (n64 > max_n)
            {
                break;
            }

            const int n =
                static_cast<int>(n64);

            if (spf[n] == n)
            {
                continue;
            }

            const Factorization f =
                factorize(n, spf);

            const vector<int> divisors =
                all_divisors(f);

            if (has_proper_one_residue_divisor(
                    n,
                    j,
                    divisors))
            {
                continue;
            }

            vector<int> exponents;

            for (size_t i = 0;
                 i < f.primes.size();
                 ++i)
            {
                const int residue =
                    f.primes[i] % j;

                const int exponent =
                    log_table[residue];

                if (exponent <= 0)
                {
                    cout << "INVALID_EXPONENT"
                         << " J=" << j
                         << " P=" << p
                         << "\n";

                    return;
                }

                for (int e = 0;
                     e < f.multiplicities[i];
                     ++e)
                {
                    exponents.push_back(exponent);
                }
            }

            int complete_sum = 0;

            for (int x : exponents)
            {
                complete_sum += x;
                complete_sum %= modulus;
            }

            if (complete_sum != 0)
            {
                cout << "ZERO_SUM_FAILURE"
                     << " J=" << j
                     << " P=" << p
                     << "\n";

                return;
            }

            if (has_proper_zero_sum_subset(
                    exponents,
                    modulus))
            {
                cout << "MINIMALITY_FAILURE"
                     << " J=" << j
                     << " P=" << p
                     << "\n";

                return;
            }

            const CanonicalSequence canonical =
                canonicalize(
                    exponents,
                    modulus
                );

            if (canonical.index <= 0)
            {
                cout << "INDEX_FAILURE"
                     << " J=" << j
                     << " P=" << p
                     << "\n";

                return;
            }

            const int omega =
                static_cast<int>(
                    canonical.values.size()
                );

            ++j_failures;
            ++total_failures;

            total_omega += omega;
            total_index += canonical.index;

            if (canonical.index == 1)
            {
                ++index_one;
            }
            else if (canonical.index == 2)
            {
                ++index_two;
                ++j_index_two;

                /*
                 * For an index-2 sequence:
                 *
                 *   sum(a_i) = 2m
                 *
                 * and minimality implies that no proper
                 * subset can sum to m.
                 */
                const bool subset_m =
                    has_subset_sum(
                        canonical.values,
                        modulus
                    );

                if (subset_m)
                {
                    ++index_two_subset_m;

                    if (omega >= 3)
                    {
                        ++index_two_subset_m_omega;
                    }

                    /*
                     * This would be surprising and is worth
                     * stopping on immediately.
                     */
                    cout << "INDEX2_SUBSET_M_FOUND"
                         << " J=" << j
                         << " P=" << p
                         << " N=" << n
                         << " OMEGA=" << omega
                         << " VALUES=";

                    print_sequence(
                        canonical.values
                    );

                    cout << "\n";
                }
                else
                {
                    ++index_two_no_subset_m;
                }

                if (j_index_two <= 8)
                {
                    cout << "INDEX2_CANONICAL"
                         << " J=" << j
                         << " P=" << p
                         << " N=" << n
                         << " MODULUS=" << modulus
                         << " OMEGA=" << omega
                         << " DISTINCT="
                         << distinct_count(
                                canonical.values
                            )
                         << " VALUES=";

                    print_sequence(
                        canonical.values
                    );

                    cout << " SUM=";

                    int sum = 0;

                    for (int x : canonical.values)
                    {
                        sum += x;
                    }

                    cout << sum
                         << " TARGET="
                         << modulus
                         << "\n";
                }
            }
            else if (canonical.index == 3)
            {
                ++index_three;
            }
            else
            {
                ++index_four_plus;
            }

            /*
             * Strong structural test:
             *
             * for canonical index 2,
             * there must be no proper subset summing m.
             */
            if (canonical.index == 2 &&
                has_subset_sum(
                    canonical.values,
                    modulus))
            {
                cout << "CONTRADICTION"
                     << " J=" << j
                     << " P=" << p
                     << " N=" << n
                     << "\n";

                return;
            }
        }

        cout << "J=" << j
             << " MODULUS=" << modulus
             << " FAILURES=" << j_failures
             << " INDEX2=" << j_index_two
             << "\n";
    }

    const double global_avg_omega =
        total_failures == 0
            ? 0.0
            : static_cast<double>(
                  total_omega
              ) /
              static_cast<double>(
                  total_failures
              );

    const double global_avg_index =
        total_failures == 0
            ? 0.0
            : static_cast<double>(
                  total_index
              ) /
              static_cast<double>(
                  total_failures
              );

    cout << "\nTOTAL_FAILURES="
         << total_failures
         << "\n";

    cout << "INDEX1="
         << index_one
         << "\n";

    cout << "INDEX2="
         << index_two
         << "\n";

    cout << "INDEX3="
         << index_three
         << "\n";

    cout << "INDEX4PLUS="
         << index_four_plus
         << "\n";

    cout << "TOTAL_INDEX="
         << total_index
         << "\n";

    cout << "TOTAL_OMEGA="
         << total_omega
         << "\n";

    cout << "GLOBAL_AVG_INDEX="
         << global_avg_index
         << "\n";

    cout << "GLOBAL_AVG_OMEGA="
         << global_avg_omega
         << "\n";

    cout << "INDEX2_SUBSET_M="
         << index_two_subset_m
         << "\n";

    cout << "INDEX2_SUBSET_M_OMEGA="
         << index_two_subset_m_omega
         << "\n";

    cout << "INDEX2_NO_SUBSET_M="
         << index_two_no_subset_m
         << "\n";

    cout << "INDEX2_NO_SUBSET_M_SANITY="
         << (
             index_two_subset_m == 0
                 ? 1
                 : 0
           )
         << "\n";

    cout << "MINIMAL_ZERO_SUM_VERIFIED=1\n";

    const auto end =
        chrono::high_resolution_clock::now();

    const double elapsed_ms =
        chrono::duration<double, milli>(
            end - start
        ).count();

    cout << "ELAPSED_TIME_MS="
         << elapsed_ms
         << "\n";
}
