# Zillion Calculator

A fully-featured **arbitrary precision scientific calculator** written in C++. Supports unlimited precision integer arithmetic, mathematical functions, and an interactive REPL interface.

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![C++](https://img.shields.io/badge/C%2B%2B-17-orange.svg)
![Platform](https://img.shields.io/badge/platform-macOS%20%7C%20Linux%20%7C%20Windows-lightgrey.svg)

## Features

### 🔢 Arbitrary Precision Arithmetic
- **Unlimited digit support** - Calculate with numbers of any size
- All basic operations: `+`, `-`, `*`, `/`, `%`, `^`
- Full **negative number support**

### 📐 Mathematical Functions
| Function | Description | Example |
|----------|-------------|---------|
| `sqrt(x)` | Square root | `sqrt(144)` → `12` |
| `cbrt(x)` | Cube root | `cbrt(27)` → `3` |
| `root(x,n)` | Nth root | `root(256,4)` → `4` |
| `abs(x)` | Absolute value | `abs(-5)` → `5` |
| `fact(n)` | Factorial | `fact(100)` → 158-digit number |
| `fib(n)` | Fibonacci | `fib(100)` → `354224848179261915075` |
| `pow(x,n)` | Power | `pow(2,100)` → 31-digit number |

### 🎲 Combinatorics
| Function | Description | Example |
|----------|-------------|---------|
| `nPr(n,r)` | Permutation | `nPr(10,3)` → `720` |
| `nCr(n,r)` | Combination | `nCr(10,5)` → `252` |

### 🔬 Number Theory
| Function | Description | Example |
|----------|-------------|---------|
| `gcd(a,b)` | Greatest Common Divisor | `gcd(48,18)` → `6` |
| `lcm(a,b)` | Least Common Multiple | `lcm(12,8)` → `24` |
| `prime(n)` | Primality test | `prime(17)` → `1` (true) |
| `factors n` | Prime factorization | `factors 1000` → `2^3 × 5^3` |

### 📊 Scientific Functions
- **Trigonometric**: `sin`, `cos`, `tan`, `asin`, `acos`, `atan`
- **Hyperbolic**: `sinh`, `cosh`, `tanh`
- **Logarithmic**: `log`, `ln`, `log10`, `exp`

### 💾 Memory Functions
| Command | Description |
|---------|-------------|
| `m+` | Add last result to memory |
| `m-` | Subtract from memory |
| `mr` | Recall memory value |
| `mc` | Clear memory |

### 🔧 Special Features
- **`ans`** - Reference previous result in calculations
- **`mem`** - Use memory value in expressions
- **Parentheses** - Full support for complex expressions
- **Continuous mode** - REPL-style interface

## Installation

### Prerequisites
- C++17 compatible compiler (g++, clang++)
- Make (optional)

### Build

```bash
# Clone the repository
git clone https://github.com/yourusername/bigINT.git
cd bigINT

# Compile
make

# Or compile directly
g++ -std=c++17 -O2 -o zillion zillion_calculator_complete.cpp
```

### Run

```bash
./zillion
```

## Usage Examples

```
>> 2^100
= 1267650600228229401496703205376 (31 digits)

>> fact(50)
= 30414093201713378043612608166064768844377641568960512000000000000 (65 digits)

>> fib(100)
= 354224848179261915075 (21 digits)

>> nCr(52,5)
= 2598960

>> gcd(123456789, 987654321)
= 9

>> sqrt(10000000000000000000000000000000000000000)
= 100000000000000000000

>> 2+3*4
= 14

>> (2+3)*4
= 20

>> ans + 100
= 120

>> m+
Memory: 120

>> factors 1000000
Prime factorization of 1000000:
  2^6 × 5^6
```

## Commands Reference

Type `help` in the calculator for the full command list.

| Command | Description |
|---------|-------------|
| `help` | Show help |
| `clear` | Clear screen |
| `digits` | Show digit count of last result |
| `factors <n>` | Prime factorization |
| `rad` | Switch to radians mode |
| `deg` | Switch to degrees mode |
| `exit` | Exit calculator |

## Technical Details

### BigInt Implementation
- Digits stored in reverse order for efficient arithmetic
- Optimized multiplication using grade-school algorithm
- Binary exponentiation for power operations
- Binary search for square/nth roots

### Expression Parser
- Recursive descent parser
- Proper operator precedence
- Support for function calls with multiple arguments
- Error handling with meaningful messages

## Performance

| Operation | 100-digit numbers | 1000-digit numbers |
|-----------|-------------------|---------------------|
| Addition | < 1ms | < 1ms |
| Multiplication | < 1ms | ~10ms |
| Division | < 1ms | ~50ms |
| Power (small exp) | < 1ms | < 10ms |

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

Created with ❤️ by Hamza

---

*Zillion Calculator - Because sometimes 64 bits just isn't enough!*
