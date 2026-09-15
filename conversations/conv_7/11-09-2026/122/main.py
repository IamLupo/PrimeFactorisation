import math
import random
import time


EXPERIMENT_ID = 148

N = 10**18 + 1234567

START = 10**8

BLOCK_SIZES = [
    256,
    1024,
    4096,
    16384,
    65536,
]

REPETITIONS = 3


def product_math_prod(n, start, size):
    end = start + size

    value = math.prod(
        range(start, end)
    )

    return value % n


def product_modular_loop(n, start, size):
    end = start + size

    value = 1

    for x in range(start, end):
        value = (
            value * x
        ) % n

    return value


def product_math_prod_chunks(
    n,
    start,
    size,
    chunk_size,
):
    end = start + size

    value = 1
    position = start

    while position < end:
        chunk_end = min(
            position + chunk_size,
            end,
        )

        chunk = math.prod(
            range(
                position,
                chunk_end,
            )
        )

        value = (
            value
            * (chunk % n)
        ) % n

        position = chunk_end

    return value


def product_tree(values, n):
    """
    Pairwise product tree.

    Every level reduces modulo n, preventing the integers from
    becoming enormous.
    """

    current = list(values)

    while len(current) > 1:
        next_level = []

        i = 0

        while i + 1 < len(current):
            next_level.append(
                (
                    current[i]
                    * current[i + 1]
                ) % n
            )

            i += 2

        if i < len(current):
            next_level.append(
                current[i] % n
            )

        current = next_level

    if current:
        return current[0] % n

    return 1


def product_tree_method(
    n,
    start,
    size,
):
    values = range(
        start,
        start + size,
    )

    return product_tree(
        values,
        n,
    )


def benchmark(
    name,
    function,
    *args,
):
    times = []
    result = None

    for _ in range(REPETITIONS):
        begin = time.perf_counter()

        result = function(
            *args
        )

        elapsed = (
            time.perf_counter()
            - begin
        )

        times.append(
            elapsed
        )

    return {
        "name": name,
        "result": result,
        "average": sum(times) / len(times),
        "minimum": min(times),
        "maximum": max(times),
    }


def main():
    print(
        f"START EXPERIMENT {EXPERIMENT_ID}"
    )
    print()

    print(
        f"N={N}"
    )
    print(
        f"START={START}"
    )
    print(
        f"repetitions={REPETITIONS}"
    )
    print()

    random.seed(
        EXPERIMENT_ID
    )

    # Small fixed test first.
    test_size = 4096

    expected = product_math_prod(
        N,
        START,
        test_size,
    )

    print("CORRECTNESS")
    print()

    methods = [
        (
            "math_prod",
            product_math_prod,
        ),
        (
            "modular_loop",
            product_modular_loop,
        ),
        (
            "tree",
            product_tree_method,
        ),
        (
            "chunks_256",
            lambda n, s, z:
            product_math_prod_chunks(
                n,
                s,
                z,
                256,
            ),
        ),
        (
            "chunks_1024",
            lambda n, s, z:
            product_math_prod_chunks(
                n,
                s,
                z,
                1024,
            ),
        ),
        (
            "chunks_4096",
            lambda n, s, z:
            product_math_prod_chunks(
                n,
                s,
                z,
                4096,
            ),
        ),
    ]

    for name, function in methods:
        result = function(
            N,
            START,
            test_size,
        )

        print(
            f"{name}: "
            f"{result == expected}"
        )

    print()

    print("BENCHMARK")
    print()

    print(
        "method\t"
        "size\t"
        "average\t"
        "minimum\t"
        "maximum"
    )

    # First benchmark a moderate block.
    for size in BLOCK_SIZES:
        print()

        print(
            f"SIZE={size}"
        )

        current_methods = [
            (
                "math_prod",
                product_math_prod,
            ),
            (
                "modular_loop",
                product_modular_loop,
            ),
            (
                "tree",
                product_tree_method,
            ),
            (
                "chunks_256",
                lambda n, s, z:
                product_math_prod_chunks(
                    n,
                    s,
                    z,
                    256,
                ),
            ),
            (
                "chunks_1024",
                lambda n, s, z:
                product_math_prod_chunks(
                    n,
                    s,
                    z,
                    1024,
                ),
            ),
            (
                "chunks_4096",
                lambda n, s, z:
                product_math_prod_chunks(
                    n,
                    s,
                    z,
                    4096,
                ),
            ),
        ]

        for name, function in current_methods:
            result = benchmark(
                name,
                function,
                N,
                START,
                size,
            )

            print(
                f"{name}\t"
                f"{size}\t"
                f"{result['average']:.6f}\t"
                f"{result['minimum']:.6f}\t"
                f"{result['maximum']:.6f}"
            )

    print()

    print(
        f"FINISHED EXPERIMENT {EXPERIMENT_ID}"
    )


if __name__ == "__main__":
    main()
