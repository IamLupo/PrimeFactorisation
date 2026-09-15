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
    int multiplier = 0;
};

struct GapRecord
{
    int delta = 0;
    int best_sum = 0;
    int best_divisor = 0;
    int complementary_divisor = 0;
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

GapRecord find_subset_gap(
    const vector<int>& values,
    int modulus,
    int n
);

int gcd_all(
    const vector<int>& values,
    int modulus
);

void print_sequence(
    const vector<int>& values
);

void print_factorization(
    const Factorization& f
);

void print_histogram(
    const unordered_map<int, long long>& histogram
);

void main_experiment();

int main()
{
    cout << "START EXPERIMENT 501\n";
    main_experiment();
    cout << "FINISHED EXPERIMENT 501\n";
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
        bool valid = true;

        for (int e : exponents)
        {
            const int x =
                (u * e) % modulus;

            if (x == 0)
            {
                valid = false;
                break;
            }

            values.push_back(x);
            sum += x;
        }

        if (!valid || sum % modulus != 0)
        {
            continue;
        }

        const int index =
            sum / modulus;

        sort(values.begin(), values.end());

        if (index < best.index ||
            (index == best.index &&
             values < best.values))
        {
            best.index = index;
            best.values = values;
            best.multiplier = u;
        }
    }

    return best;
}

GapRecord find_subset_gap(
    const vector<int>& values,
    int modulus,
    int n
)
{
    GapRecord result;

    vector<bool> possible(modulus, false);
    vector<int> best_subset(modulus, 0);

    possible[0] = true;

    for (int x : values)
    {
        vector<bool> next = possible;

        for (int s = 0; s < modulus; ++s)
        {
            if (!possible[s])
            {
                continue;
            }

            const int ns = s + x;

            if (ns <= modulus - 1 &&
                !next[ns])
            {
                next[ns] = true;
                best_subset[ns] = best_subset[s] + x;
            }
        }

        possible.swap(next);
    }

    /*
     * We are interested in the largest subset sum strictly
     * below modulus.
     */
    int best_sum = 0;

    for (int s = 1; s < modulus; ++s)
    {
        if (possible[s])
        {
            best_sum = s;
        }
    }

    result.best_sum = best_sum;
    result.delta = modulus - best_sum;

    /*
     * Reconstruct an actual divisor is not possible from only
     * the canonical exponent values, because the canonical
     * multiplier changes the exponent representation.
     *
     * The caller therefore uses the residue value
     * g^(m-delta) separately.
     */
    result.best_divisor = 0;
    result.complementary_divisor = 0;

    return result;
}

int gcd_all(
    const vector<int>& values,
    int modulus
)
{
    int result = modulus;

    for (int x : values)
    {
        result = gcd(result, x);
    }

    return result;
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

void print_histogram(
    const unordered_map<int, long long>& histogram
)
{
    vector<pair<int, long long>> values(
        histogram.begin(),
        histogram.end()
    );

    sort(
        values.begin(),
        values.end(),
        [](const auto& a, const auto& b)
        {
            return a.first < b.first;
        }
    );

    for (const auto& [key, value] : values)
    {
        cout << key << ":" << value << " ";
    }

    cout << "\n";
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
         << primes.size() << "\n";

    cout << "J_LIMIT="
         << J_LIMIT << "\n";

    cout << "MAX_JP_PLUS_1="
         << max_n << "\n";

    cout << "SPF_READY=1\n";

    long long total_index2 = 0;
    long long total_index2_gap = 0;

    unordered_map<int, long long> gap_hist;
    unordered_map<int, long long> omega_hist;

    int global_min_delta = 1000000000;
    int global_max_delta = 0;

    for (int j = 2; j <= J_LIMIT; ++j)
    {
        if (!is_prime(j))
        {
            continue;
        }

        const int modulus = j - 1;

        if (modulus <= 1)
        {
            continue;
        }

        const int primitive_root =
            primitive_root_prime(j);

        if (primitive_root < 0)
        {
            cout << "PRIMITIVE_ROOT_FAILURE J="
                 << j
                 << "\n";

            return;
        }

        const vector<int> logs =
            build_discrete_log_table(
                primitive_root,
                j
            );

        long long j_index2 = 0;

        unordered_map<int, long long>
            j_gap_hist;

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

                const int e =
                    logs[residue];

                for (int k = 0;
                     k < f.multiplicities[i];
                     ++k)
                {
                    exponents.push_back(e);
                }
            }

            if (has_proper_zero_sum_subset(
                    exponents,
                    modulus))
            {
                cout << "MINIMALITY_FAILURE\n";
                return;
            }

            const CanonicalSequence canonical =
                canonicalize(
                    exponents,
                    modulus
                );

            if (canonical.index != 2)
            {
                continue;
            }

            ++j_index2;
            ++total_index2;

            const GapRecord gap =
                find_subset_gap(
                    canonical.values,
                    modulus,
                    n
                );

            ++total_index2_gap;
            ++gap_hist[gap.delta];
            ++j_gap_hist[gap.delta];

            const int omega =
                static_cast<int>(
                    canonical.values.size()
                );

            ++omega_hist[omega];

            global_min_delta =
                min(
                    global_min_delta,
                    gap.delta
                );

            global_max_delta =
                max(
                    global_max_delta,
                    gap.delta
                );

            if (j_index2 <= 8)
            {
                const int near_residue =
                    mod_pow(
                        primitive_root,
                        modulus - gap.delta,
                        j
                    );

                cout << "INDEX2_GAP"
                     << " J=" << j
                     << " P=" << p
                     << " N=" << n
                     << " MODULUS="
                     << modulus
                     << " OMEGA="
                     << omega
                     << " DELTA="
                     << gap.delta
                     << " BEST_SUM="
                     << gap.best_sum
                     << " NEAR_RESIDUE="
                     << near_residue
                     << " VALUES=";

                print_sequence(
                    canonical.values
                );

                cout << "\n";
            }
        }

        cout << "J=" << j
             << " MODULUS=" << modulus
             << " INDEX2="
             << j_index2
             << " GAP_HIST=";

        print_histogram(j_gap_hist);
    }

    cout << "\nTOTAL_INDEX2="
         << total_index2
         << "\n";

    cout << "TOTAL_INDEX2_GAP="
         << total_index2_gap
         << "\n";

    cout << "MIN_DELTA="
         << global_min_delta
         << "\n";

    cout << "MAX_DELTA="
         << global_max_delta
         << "\n";

    cout << "GAP_HIST=";
    print_histogram(gap_hist);

    cout << "INDEX2_OMEGA_HIST=";
    print_histogram(omega_hist);

    cout << "GAP_SANITY="
         << (
             total_index2 ==
             total_index2_gap
                 ? 1
                 : 0
           )
         << "\n";

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
