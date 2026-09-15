#include <cstdint>
#include <iostream>
#include <chrono>
#include <numeric>

using u64 = std::uint64_t;
using u128 = unsigned __int128;


u64 mul_mod(u64 a, u64 b, u64 n)
{
    return static_cast<u64>(
        (static_cast<u128>(a) * b) % n
    );
}


u64 range_product_mod(
    u64 start,
    u64 end,
    u64 n
)
{
    u64 result = 1;

    for (u64 x = start; x <= end; ++x)
    {
        result = mul_mod(
            result,
            x,
            n
        );
    }

    return result;
}


u64 range_product_blocked(
    u64 start,
    u64 end,
    u64 n,
    u64 block_size
)
{
    u64 result = 1;

    u64 position = start;

    while (position <= end)
    {
        u64 block_end =
            std::min(
                position + block_size - 1,
                end
            );

        u64 block_product = 1;

        for (u64 x = position;
             x <= block_end;
             ++x)
        {
            block_product = mul_mod(
                block_product,
                x,
                n
            );
        }

        result = mul_mod(
            result,
            block_product,
            n
        );

        position = block_end + 1;
    }

    return result;
}


bool benchmark(
    const char *name,
    u64 start,
    u64 end,
    u64 n,
    u64 block_size
)
{
    auto begin =
        std::chrono::steady_clock::now();

    u64 result;

    if (block_size == 0)
    {
        result = range_product_mod(
            start,
            end,
            n
        );
    }
    else
    {
        result = range_product_blocked(
            start,
            end,
            n,
            block_size
        );
    }

    auto finish =
        std::chrono::steady_clock::now();

    double seconds =
        std::chrono::duration<double>(
            finish - begin
        ).count();

    std::cout
        << name
        << " result="
        << result
        << " seconds="
        << seconds
        << "\n";

    return true;
}


int main()
{
    std::cout
        << "START EXPERIMENT 150\n\n";

    // Same scale as the successful 10^18 Python case.
    const u64 n =
        1000000005045401363ULL;

    const u64 start =
        948054019ULL;

    const u64 end =
        1000000002ULL;

    std::cout
        << "N=" << n << "\n"
        << "start=" << start << "\n"
        << "end=" << end << "\n"
        << "terms="
        << (end - start + 1)
        << "\n\n";

    // Correctness reference.
    u64 reference =
        range_product_blocked(
            start,
            start + 4095,
            n,
            256
        );

    u64 direct =
        range_product_mod(
            start,
            start + 4095,
            n
        );

    std::cout
        << "CORRECTNESS="
        << (reference == direct)
        << "\n\n";

    benchmark(
        "direct_mod_loop",
        start,
        end,
        n,
        0
    );

    benchmark(
        "block_256",
        start,
        end,
        n,
        256
    );

    benchmark(
        "block_1024",
        start,
        end,
        n,
        1024
    );

    benchmark(
        "block_4096",
        start,
        end,
        n,
        4096
    );

    benchmark(
        "block_16384",
        start,
        end,
        n,
        16384
    );

    std::cout
        << "\nFINISHED EXPERIMENT 150\n";

    return 0;
}
