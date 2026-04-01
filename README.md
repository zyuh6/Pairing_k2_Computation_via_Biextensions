# Computing Pairings on Elliptic Curves with Embedding Degree Two via Biextensions

**Paper:** [Computing Pairings on Elliptic Curves with Embedding Degree Two via Biextensions](https://eprint.iacr.org/2025/1652)

This repository contains the implementation accompanying the research paper *"Computing Pairings on Elliptic Curves with Embedding Degree Two via Biextensions"*. It provides the source code to implement, verify, and benchmark the performance (in CPU clock cycles) of six distinct pairing computations on supersingular (SS) and non-supersingular (NSS) elliptic curves with an embedding degree of 2 over 1536-bit fields.

## Implemented Pairings

| ID | Pairing Type | Curve Profile | Algorithm Used |
|:---:|:---|:---|:---|
| 1 | Tate Pairing | `SSmont1536`| Double-and-add ladder with line function |
| 2 | Variant of Tate Pairing | `NSS1536` | Double-and-add ladder with line function |
| 3 | Omega Pairing | `NSS1536` | Double-and-add ladder with line function |
| 4 | Tate Pairing | `SS1536` | Miller's algorithm |
| 5 | Variant of Tate Pairing | `NSS1536(a=-3)` | Miller's algorithm |
| 6 | Omega Pairing | `NSS1536(a=-3)` | Miller's algorithm |

## Project Structure

The source code is organized into two separate directories based on curve parameters. Each requires an independent build process:

* `relic-ss1536_and_nss1536_aeq1/`: Pairings 1-4 (SS curves and NSS curves with `a=1`).
* `relic-nss1536_aeqn3/`: Pairings 5-6 (NSS curves with `a=-3`).

## Dependencies

This project relies on a modified version of the [RELIC cryptography toolkit](https://github.com/relic-toolkit/relic). Ensure the following build tools are installed on your system: CMake, Make, GMP, and a compatible C compiler (e.g., gcc or clang).

## Building and Testing

The following commands use standard indentation to represent terminal inputs.

### Part A: Pairings 1-4 (`a=1` and SS)

Navigate to the directory:
    
    cd relic-ss1536_and_nss1536_aeq1

Option 1: Build for SS1536 (Pairings 1 & 4)
    
    mkdir build_ss && cd build_ss
    ../preset/gmp-pbc-ss1536.sh ../
    make test_1536
    ./bin/test_1536
    make bench_1536
    ./bin/bench_1536

Option 2: Build for NSS1536 (Pairings 2 & 3)
*(Run this from the root of `relic-ss1536_and_nss1536_aeq1`)*
    
    mkdir build_nss && cd build_nss
    ../preset/gmp-pbc-nss1536.sh ../
    make test_1536
    ./bin/test_1536
    make bench_1536
    ./bin/bench_1536

### Part B: Pairings 5-6 (`a=-3`)

Navigate to the directory:
    
    cd relic-nss1536_aeqn3

Build and Evaluate:
    
    mkdir build && cd build
    ../preset/gmp-pbc-nss1536.sh ../
    make test_1536
    ./bin/test_1536
    make bench_1536
    ./bin/bench_1536

## Results

* **Correctness Verification (`test_1536`)**: Outputs cryptographic assertions to verify the correct computation of the pairings over the respective algebraic structures.
* **Performance Benchmarks (`bench_1536`)**: Measures and outputs the execution time in CPU clock cycles.

## License & Copyright

The source code in this repository is licensed under the **[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0)**, see the LICENSE file for details.

This implementation is built upon **[RELIC cryptography toolkit](https://github.com/relic-toolkit/relic)**, which is dual-licensed under LGPL 2.1+ and Apache 2.0.

The accompanying research paper and its associated figures/data are licensed under CC BY 4.0.

