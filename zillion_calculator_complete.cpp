/**
 * ╔═══════════════════════════════════════════════════════════════════════════╗
 * ║              ZILLION SCIENTIFIC CALCULATOR - BigInt Library               ║
 * ║                                                                           ║
 * ║  A fully-fledged scientific calculator with arbitrary precision:         ║
 * ║  • Continuous calculation mode (like a real calculator)                  ║
 * ║  • All arithmetic operators with unlimited precision                     ║
 * ║  • Scientific functions (sqrt, log, sin, cos, tan, etc.)                ║
 * ║  • Mathematical functions (factorial, fibonacci, permutation, etc.)     ║
 * ║  • Memory functions (M+, M-, MR, MC)                                    ║
 * ║  • Expression evaluation with function calls                            ║
 * ║  • Number theory (GCD, LCM, primality, prime factorization)            ║
 * ║  • Constants (pi, e, golden ratio)                                      ║
 * ╚═══════════════════════════════════════════════════════════════════════════╝
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

// ═══════════════════════════════════════════════════════════════════════════
//  CROSS-PLATFORM COMPATIBILITY
// ═══════════════════════════════════════════════════════════════════════════

#ifdef _WIN32
#define CLEAR_SCREEN "cls"
#include <windows.h>
void sleepMs(int ms) { Sleep(ms); }
#else
#define CLEAR_SCREEN "clear"
#include <unistd.h>
void sleepMs(int ms) { usleep(ms * 1000); }
#endif

// ═══════════════════════════════════════════════════════════════════════════
//  BIGINT CLASS - ARBITRARY PRECISION INTEGERS
// ═══════════════════════════════════════════════════════════════════════════

class BigInt {
private:
  string digits;
  bool is_negative;

public:
  BigInt(long long n = 0);
  BigInt(const string &s);

  friend void divide_by_2(BigInt &a);
  friend bool Null(const BigInt &a);
  friend int Length(const BigInt &a);
  friend BigInt Abs(const BigInt &a);
  int operator[](int index) const;

  BigInt &operator=(const BigInt &other);
  BigInt operator-() const;
  BigInt &operator++();
  BigInt operator++(int);
  BigInt &operator--();
  BigInt operator--(int);

  friend bool operator==(const BigInt &a, const BigInt &b);
  friend bool operator!=(const BigInt &a, const BigInt &b);
  friend bool operator<(const BigInt &a, const BigInt &b);
  friend bool operator<=(const BigInt &a, const BigInt &b);
  friend bool operator>(const BigInt &a, const BigInt &b);
  friend bool operator>=(const BigInt &a, const BigInt &b);

  friend BigInt operator+(const BigInt &a, const BigInt &b);
  friend BigInt operator-(const BigInt &a, const BigInt &b);
  friend BigInt operator*(const BigInt &a, const BigInt &b);
  friend BigInt operator/(const BigInt &a, const BigInt &b);
  friend BigInt operator%(const BigInt &a, const BigInt &b);
  friend BigInt power(const BigInt &base, const BigInt &exp);

  friend BigInt &operator+=(BigInt &a, const BigInt &b);
  friend BigInt &operator-=(BigInt &a, const BigInt &b);
  friend BigInt &operator*=(BigInt &a, const BigInt &b);
  friend BigInt &operator/=(BigInt &a, const BigInt &b);
  friend BigInt &operator%=(BigInt &a, const BigInt &b);

  friend istream &operator>>(istream &in, BigInt &a);
  friend ostream &operator<<(ostream &out, const BigInt &a);

  bool isNegative() const { return is_negative; }
  bool isZero() const { return Null(*this); }
  string toString() const;
  int digitCount() const { return digits.size(); }
  double toDouble() const;

private:
  friend BigInt &addMagnitudes(BigInt &a, const BigInt &b);
  friend BigInt &subtractMagnitudes(BigInt &a, const BigInt &b);
  friend BigInt &multiplyMagnitudes(BigInt &a, const BigInt &b);
  friend BigInt &divideMagnitudes(BigInt &a, const BigInt &b);
  friend BigInt &incrementMagnitude(BigInt &a);
  friend BigInt &decrementMagnitude(BigInt &a);
  void removeLeadingZeros();
};

// ─────────────────────────── BigInt Implementation ───────────────────────────

BigInt::BigInt(long long n) : is_negative(false) {
  if (n < 0) {
    is_negative = true;
    n = -n;
  }
  do {
    digits.push_back(n % 10);
    n /= 10;
  } while (n);
}

BigInt::BigInt(const string &s) : is_negative(false) {
  digits = "";
  int startIdx = 0;
  if (s.empty()) {
    digits.push_back(0);
    return;
  }
  if (s[0] == '-') {
    is_negative = true;
    startIdx = 1;
  } else if (s[0] == '+') {
    startIdx = 1;
  }
  for (int i = s.length() - 1; i >= startIdx; i--) {
    if (s[i] < '0' || s[i] > '9')
      throw invalid_argument("Invalid number: " + s);
    digits.push_back(s[i] - '0');
  }
  removeLeadingZeros();
  if (Null(*this))
    is_negative = false;
}

void BigInt::removeLeadingZeros() {
  while (digits.size() > 1 && digits.back() == 0)
    digits.pop_back();
}

bool Null(const BigInt &a) { return a.digits.size() == 1 && a.digits[0] == 0; }
int Length(const BigInt &a) { return a.digits.size(); }
BigInt Abs(const BigInt &a) {
  BigInt t = a;
  t.is_negative = false;
  return t;
}

int BigInt::operator[](int index) const {
  if (index < 0 || index >= (int)digits.size())
    throw out_of_range("Index out of bounds");
  return digits[index];
}

void divide_by_2(BigInt &a) {
  int carry = 0;
  for (int i = a.digits.size() - 1; i >= 0; i--) {
    int digit = (a.digits[i] >> 1) + carry;
    carry = (a.digits[i] & 1) * 5;
    a.digits[i] = digit;
  }
  a.removeLeadingZeros();
}

BigInt &BigInt::operator=(const BigInt &o) {
  if (this != &o) {
    digits = o.digits;
    is_negative = o.is_negative;
  }
  return *this;
}

BigInt BigInt::operator-() const {
  BigInt t = *this;
  if (!Null(t))
    t.is_negative = !t.is_negative;
  return t;
}

BigInt &incrementMagnitude(BigInt &a) {
  int i = 0, n = a.digits.size();
  while (i < n && a.digits[i] == 9)
    a.digits[i++] = 0;
  if (i == n)
    a.digits.push_back(1);
  else
    a.digits[i]++;
  return a;
}

BigInt &decrementMagnitude(BigInt &a) {
  int i = 0, n = a.digits.size();
  while (i < n && a.digits[i] == 0)
    a.digits[i++] = 9;
  a.digits[i]--;
  a.removeLeadingZeros();
  return a;
}

BigInt &BigInt::operator++() {
  if (is_negative) {
    decrementMagnitude(*this);
    if (Null(*this))
      is_negative = false;
  } else
    incrementMagnitude(*this);
  return *this;
}

BigInt BigInt::operator++(int) {
  BigInt t = *this;
  ++(*this);
  return t;
}

BigInt &BigInt::operator--() {
  if (Null(*this)) {
    incrementMagnitude(*this);
    is_negative = true;
  } else if (is_negative)
    incrementMagnitude(*this);
  else
    decrementMagnitude(*this);
  return *this;
}

BigInt BigInt::operator--(int) {
  BigInt t = *this;
  --(*this);
  return t;
}

bool operator==(const BigInt &a, const BigInt &b) {
  return a.digits == b.digits && a.is_negative == b.is_negative;
}
bool operator!=(const BigInt &a, const BigInt &b) { return !(a == b); }

bool operator<(const BigInt &a, const BigInt &b) {
  if (a.is_negative && !b.is_negative)
    return true;
  if (!a.is_negative && b.is_negative)
    return false;
  int n = Length(a), m = Length(b);
  if (n != m)
    return a.is_negative ? (n > m) : (n < m);
  for (int i = n - 1; i >= 0; i--) {
    if (a.digits[i] != b.digits[i]) {
      bool less = a.digits[i] < b.digits[i];
      return a.is_negative ? !less : less;
    }
  }
  return false;
}

bool operator<=(const BigInt &a, const BigInt &b) { return !(a > b); }
bool operator>(const BigInt &a, const BigInt &b) { return b < a; }
bool operator>=(const BigInt &a, const BigInt &b) { return !(a < b); }

BigInt &addMagnitudes(BigInt &a, const BigInt &b) {
  int n = Length(a), m = Length(b);
  if (n < m) {
    a.digits.append(m - n, 0);
    n = m;
  }
  int carry = 0;
  for (int i = 0; i < n; i++) {
    int sum = a.digits[i] + (i < m ? b.digits[i] : 0) + carry;
    a.digits[i] = sum % 10;
    carry = sum / 10;
  }
  if (carry)
    a.digits.push_back(carry);
  return a;
}

BigInt &subtractMagnitudes(BigInt &a, const BigInt &b) {
  int n = Length(a), m = Length(b), borrow = 0;
  for (int i = 0; i < n; i++) {
    int diff = a.digits[i] - (i < m ? b.digits[i] : 0) + borrow;
    if (diff < 0) {
      diff += 10;
      borrow = -1;
    } else
      borrow = 0;
    a.digits[i] = diff;
  }
  a.removeLeadingZeros();
  return a;
}

BigInt &multiplyMagnitudes(BigInt &a, const BigInt &b) {
  if (Null(a) || Null(b)) {
    a = BigInt(0LL);
    return a;
  }
  int n = a.digits.size(), m = b.digits.size();
  vector<int> result(n + m, 0);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      result[i + j] += a.digits[i] * b.digits[j];
  int carry = 0;
  for (int i = 0; i < n + m; i++) {
    int sum = result[i] + carry;
    result[i] = sum % 10;
    carry = sum / 10;
  }
  a.digits.resize(n + m);
  for (int i = 0; i < n + m; i++)
    a.digits[i] = result[i];
  a.removeLeadingZeros();
  return a;
}

BigInt &divideMagnitudes(BigInt &a, const BigInt &b) {
  if (Null(b))
    throw invalid_argument("Division by zero");
  BigInt absA = Abs(a), absB = Abs(b);
  if (absA < absB) {
    a = BigInt(0LL);
    return a;
  }
  if (absA == absB) {
    a = BigInt(1LL);
    return a;
  }
  BigInt t(0LL);
  vector<int> quotient;
  int n = Length(a), i = n - 1;
  while (t * 10 + a.digits[i] < absB && i > 0) {
    t = t * 10 + a.digits[i];
    i--;
  }
  for (; i >= 0; i--) {
    t = t * 10 + a.digits[i];
    int q;
    for (q = 9; q * absB > t; q--)
      ;
    t -= q * absB;
    quotient.push_back(q);
  }
  int len = quotient.size();
  a.digits.resize(len);
  for (int j = 0; j < len; j++)
    a.digits[j] = quotient[len - j - 1];
  a.removeLeadingZeros();
  return a;
}

BigInt &operator+=(BigInt &a, const BigInt &b) {
  BigInt absA = Abs(a), absB = Abs(b);
  if (a.is_negative == b.is_negative) {
    bool neg = a.is_negative;
    addMagnitudes(a, b);
    a.is_negative = neg;
  } else if (absA > absB) {
    bool neg = a.is_negative;
    BigInt t = absA;
    subtractMagnitudes(t, absB);
    a = t;
    a.is_negative = neg;
  } else if (absA < absB) {
    bool neg = b.is_negative;
    BigInt t = absB;
    subtractMagnitudes(t, absA);
    a = t;
    a.is_negative = neg;
  } else
    a = BigInt(0LL);
  if (Null(a))
    a.is_negative = false;
  return a;
}

BigInt operator+(const BigInt &a, const BigInt &b) {
  BigInt r = a;
  r += b;
  return r;
}
BigInt &operator-=(BigInt &a, const BigInt &b) {
  a += (-b);
  return a;
}
BigInt operator-(const BigInt &a, const BigInt &b) {
  BigInt r = a;
  r -= b;
  return r;
}

BigInt &operator*=(BigInt &a, const BigInt &b) {
  bool neg = a.is_negative != b.is_negative;
  BigInt absA = Abs(a), absB = Abs(b);
  multiplyMagnitudes(absA, absB);
  a = absA;
  a.is_negative = neg && !Null(a);
  return a;
}

BigInt operator*(const BigInt &a, const BigInt &b) {
  BigInt r = a;
  r *= b;
  return r;
}

BigInt &operator/=(BigInt &a, const BigInt &b) {
  bool neg = a.is_negative != b.is_negative;
  BigInt absA = Abs(a), absB = Abs(b);
  divideMagnitudes(absA, absB);
  a = absA;
  a.is_negative = neg && !Null(a);
  return a;
}

BigInt operator/(const BigInt &a, const BigInt &b) {
  BigInt r = a;
  r /= b;
  return r;
}

BigInt &operator%=(BigInt &a, const BigInt &b) {
  if (Null(b))
    throw invalid_argument("Modulo by zero");
  BigInt absA = Abs(a), absB = Abs(b);
  if (absA < absB)
    return a;
  if (absA == absB) {
    a = BigInt(0LL);
    return a;
  }
  BigInt t(0LL);
  int n = Length(a), i = n - 1;
  while (t * 10 + a.digits[i] < absB && i > 0)
    t = t * 10 + a.digits[i--];
  for (; i >= 0; i--) {
    t = t * 10 + a.digits[i];
    int q;
    for (q = 9; q * absB > t; q--)
      ;
    t -= q * absB;
  }
  bool neg = a.is_negative;
  a = t;
  a.is_negative = neg && !Null(a);
  return a;
}

BigInt operator%(const BigInt &a, const BigInt &b) {
  BigInt r = a;
  r %= b;
  return r;
}

BigInt power(const BigInt &base, const BigInt &exp) {
  if (exp.is_negative)
    return BigInt(0LL);
  BigInt b = base, e = exp, result(1LL);
  while (!Null(e)) {
    if (e.digits[0] & 1)
      result *= b;
    b *= b;
    divide_by_2(e);
  }
  return result;
}

istream &operator>>(istream &in, BigInt &a) {
  string s;
  in >> s;
  a = BigInt(s);
  return in;
}

ostream &operator<<(ostream &out, const BigInt &a) {
  if (a.is_negative)
    out << "-";
  for (int i = a.digits.size() - 1; i >= 0; i--)
    out << (int)a.digits[i];
  return out;
}

string BigInt::toString() const {
  string r;
  if (is_negative)
    r = "-";
  for (int i = digits.size() - 1; i >= 0; i--)
    r += ('0' + digits[i]);
  return r;
}

double BigInt::toDouble() const {
  double result = 0, multiplier = 1;
  for (size_t i = 0; i < digits.size(); i++) {
    result += digits[i] * multiplier;
    multiplier *= 10;
  }
  return is_negative ? -result : result;
}

// ═══════════════════════════════════════════════════════════════════════════
//  MATHEMATICAL FUNCTIONS FOR BIGINT
// ═══════════════════════════════════════════════════════════════════════════

BigInt bigSqrt(BigInt a) {
  if (a.isNegative())
    throw invalid_argument("Square root of negative number");
  if (Null(a) || a == BigInt(1LL))
    return a;
  BigInt left(1LL), right(a), result(0LL), mid;
  while (left <= right) {
    mid = (left + right) / 2;
    BigInt sq = mid * mid;
    if (sq == a)
      return mid;
    else if (sq < a) {
      result = mid;
      left = mid + 1;
    } else
      right = mid - 1;
  }
  return result;
}

BigInt nthRoot(BigInt a, BigInt n) {
  if (a.isNegative() && n[0] % 2 == 0)
    throw invalid_argument("Even root of negative");
  if (Null(a))
    return BigInt(0LL);
  if (a == BigInt(1LL) || n == BigInt(1LL))
    return a;
  BigInt low(1LL), high = a, mid, result;
  while (low <= high) {
    mid = (low + high) / 2;
    BigInt p = power(mid, n);
    if (p == a)
      return mid;
    else if (p < a) {
      result = mid;
      low = mid + 1;
    } else
      high = mid - 1;
  }
  return result;
}

BigInt gcd(BigInt a, BigInt b) {
  a = Abs(a);
  b = Abs(b);
  while (!Null(b)) {
    BigInt t = b;
    b = a % b;
    a = t;
  }
  return a;
}

BigInt lcm(const BigInt &a, const BigInt &b) {
  return (Abs(a) / gcd(a, b)) * Abs(b);
}

BigInt factorial(BigInt n) {
  if (n.isNegative())
    throw invalid_argument("Factorial of negative");
  BigInt result(1LL);
  for (BigInt i = 2LL; i <= n; ++i)
    result *= i;
  return result;
}

BigInt fibonacci(int n) {
  if (n < 0)
    throw invalid_argument("Fibonacci of negative index");
  if (n == 0)
    return BigInt(0LL);
  if (n == 1)
    return BigInt(1LL);
  BigInt a(0LL), b(1LL), c;
  for (int i = 2; i <= n; i++) {
    c = a + b;
    a = b;
    b = c;
  }
  return b;
}

BigInt permutation(BigInt n, BigInt r) {
  if (r > n)
    return BigInt(0LL);
  BigInt result(1LL);
  for (BigInt i = 0LL; i < r; ++i)
    result *= (n - i);
  return result;
}

BigInt combination(BigInt n, BigInt r) {
  if (r > n)
    return BigInt(0LL);
  if (r > n - r)
    r = n - r;
  BigInt result(1LL);
  for (BigInt i = 0LL; i < r; ++i) {
    result *= (n - i);
    result /= (i + 1);
  }
  return result;
}

bool isPrime(BigInt n) {
  if (n < BigInt(2LL))
    return false;
  if (n == BigInt(2LL))
    return true;
  if (n % 2 == BigInt(0LL))
    return false;
  BigInt sqrtN = bigSqrt(n);
  for (BigInt i = 3LL; i <= sqrtN; i += 2LL)
    if (n % i == BigInt(0LL))
      return false;
  return true;
}

vector<pair<BigInt, int>> primeFactorize(BigInt n) {
  vector<pair<BigInt, int>> factors;
  n = Abs(n);
  int count = 0;
  while (n % 2 == BigInt(0LL)) {
    count++;
    n /= 2;
  }
  if (count > 0)
    factors.push_back({BigInt(2LL), count});
  for (BigInt i = 3LL; i * i <= n; i += 2LL) {
    count = 0;
    while (n % i == BigInt(0LL)) {
      count++;
      n /= i;
    }
    if (count > 0)
      factors.push_back({i, count});
  }
  if (n > BigInt(1LL))
    factors.push_back({n, 1});
  return factors;
}

// ═══════════════════════════════════════════════════════════════════════════
//  SCIENTIFIC FUNCTIONS (FLOATING POINT)
// ═══════════════════════════════════════════════════════════════════════════

const double PI = 3.14159265358979323846;
const double E = 2.71828182845904523536;
const double PHI = 1.61803398874989484820; // Golden ratio

double toRadians(double degrees) { return degrees * PI / 180.0; }
double toDegrees(double radians) { return radians * 180.0 / PI; }

// ═══════════════════════════════════════════════════════════════════════════
//  CALCULATOR STATE
// ═══════════════════════════════════════════════════════════════════════════

struct CalculatorState {
  BigInt lastResult;
  BigInt memory;
  bool hasLastResult;
  bool useRadians;
  int precision;

  CalculatorState()
      : lastResult(0LL), memory(0LL), hasLastResult(false), useRadians(true),
        precision(10) {}
};

CalculatorState calcState;

// ═══════════════════════════════════════════════════════════════════════════
//  EXPRESSION PARSER WITH FUNCTION SUPPORT
// ═══════════════════════════════════════════════════════════════════════════

string trim(const string &s) {
  size_t start = s.find_first_not_of(" \t\n\r");
  size_t end = s.find_last_not_of(" \t\n\r");
  return (start == string::npos) ? "" : s.substr(start, end - start + 1);
}

string toLower(const string &s) {
  string result = s;
  transform(result.begin(), result.end(), result.begin(), ::tolower);
  return result;
}

bool isOperator(char c) {
  return c == '+' || c == '-' || c == '*' || c == '/' || c == '%' || c == '^';
}

int precedence(char op) {
  if (op == '^')
    return 3;
  if (op == '*' || op == '/' || op == '%')
    return 2;
  if (op == '+' || op == '-')
    return 1;
  return 0;
}

bool isFunction(const string &name) {
  static const vector<string> funcs = {
      "sqrt", "cbrt",  "abs", "fact", "fib",   "gcd",  "lcm",
      "sin",  "cos",   "tan", "asin", "acos",  "atan", "sinh",
      "cosh", "tanh",  "log", "ln",   "log10", "exp",  "floor",
      "ceil", "round", "npr", "ncr",  "prime", "root"};
  string lower = toLower(name);
  for (const auto &f : funcs)
    if (f == lower)
      return true;
  return false;
}

class ExpressionParser {
private:
  string expr;
  size_t pos;

  char peek() { return pos < expr.size() ? expr[pos] : '\0'; }
  char get() { return pos < expr.size() ? expr[pos++] : '\0'; }
  void skip() {
    while (pos < expr.size() && isspace(expr[pos]))
      pos++;
  }

  BigInt parseNumber() {
    skip();
    string num;
    bool neg = false;
    if (peek() == '-') {
      neg = true;
      get();
      skip();
    }
    while (isdigit(peek()))
      num += get();
    if (num.empty())
      throw invalid_argument("Expected number");
    return neg ? -BigInt(num) : BigInt(num);
  }

  string parseIdentifier() {
    skip();
    string id;
    while (isalpha(peek()) || peek() == '_')
      id += get();
    return id;
  }

  BigInt parsePrimary() {
    skip();

    // Handle 'ans' for last result
    if (toLower(expr.substr(pos, 3)) == "ans") {
      pos += 3;
      if (!calcState.hasLastResult)
        throw invalid_argument("No previous result");
      return calcState.lastResult;
    }

    // Handle 'mem' for memory
    if (toLower(expr.substr(pos, 3)) == "mem") {
      pos += 3;
      return calcState.memory;
    }

    // Handle constants
    if (toLower(expr.substr(pos, 2)) == "pi") {
      pos += 2;
      return BigInt(3LL); // Integer approximation
    }

    // Handle parentheses
    if (peek() == '(') {
      get();
      BigInt result = parseExpression();
      skip();
      if (peek() != ')')
        throw invalid_argument("Missing closing parenthesis");
      get();
      return result;
    }

    // Handle unary minus
    if (peek() == '-') {
      get();
      return -parsePrimary();
    }

    // Handle unary plus
    if (peek() == '+') {
      get();
      return parsePrimary();
    }

    // Handle functions
    if (isalpha(peek())) {
      string funcName = parseIdentifier();
      string lower = toLower(funcName);
      skip();

      if (peek() == '(') {
        get();
        vector<BigInt> args;

        if (peek() != ')') {
          args.push_back(parseExpression());
          skip();
          while (peek() == ',') {
            get();
            args.push_back(parseExpression());
            skip();
          }
        }

        if (peek() != ')')
          throw invalid_argument("Missing closing parenthesis");
        get();

        return evaluateFunction(lower, args);
      }

      throw invalid_argument("Unknown identifier: " + funcName);
    }

    // Handle numbers
    if (isdigit(peek())) {
      return parseNumber();
    }

    throw invalid_argument("Unexpected character: " + string(1, peek()));
  }

  BigInt evaluateFunction(const string &name, const vector<BigInt> &args) {
    if (name == "sqrt") {
      if (args.size() != 1)
        throw invalid_argument("sqrt requires 1 argument");
      return bigSqrt(args[0]);
    }
    if (name == "cbrt" || name == "root") {
      if (args.size() == 1)
        return nthRoot(args[0], BigInt(3LL));
      if (args.size() == 2)
        return nthRoot(args[0], args[1]);
      throw invalid_argument("root requires 1 or 2 arguments");
    }
    if (name == "abs") {
      if (args.size() != 1)
        throw invalid_argument("abs requires 1 argument");
      return Abs(args[0]);
    }
    if (name == "fact" || name == "factorial") {
      if (args.size() != 1)
        throw invalid_argument("fact requires 1 argument");
      return factorial(args[0]);
    }
    if (name == "fib" || name == "fibonacci") {
      if (args.size() != 1)
        throw invalid_argument("fib requires 1 argument");
      return fibonacci((int)args[0].toDouble());
    }
    if (name == "gcd") {
      if (args.size() != 2)
        throw invalid_argument("gcd requires 2 arguments");
      return gcd(args[0], args[1]);
    }
    if (name == "lcm") {
      if (args.size() != 2)
        throw invalid_argument("lcm requires 2 arguments");
      return lcm(args[0], args[1]);
    }
    if (name == "npr" || name == "perm") {
      if (args.size() != 2)
        throw invalid_argument("nPr requires 2 arguments");
      return permutation(args[0], args[1]);
    }
    if (name == "ncr" || name == "comb") {
      if (args.size() != 2)
        throw invalid_argument("nCr requires 2 arguments");
      return combination(args[0], args[1]);
    }
    if (name == "prime" || name == "isprime") {
      if (args.size() != 1)
        throw invalid_argument("prime requires 1 argument");
      return isPrime(args[0]) ? BigInt(1LL) : BigInt(0LL);
    }
    if (name == "pow" || name == "power") {
      if (args.size() != 2)
        throw invalid_argument("pow requires 2 arguments");
      return power(args[0], args[1]);
    }
    if (name == "mod") {
      if (args.size() != 2)
        throw invalid_argument("mod requires 2 arguments");
      return args[0] % args[1];
    }
    if (name == "min") {
      if (args.size() < 2)
        throw invalid_argument("min requires at least 2 arguments");
      BigInt result = args[0];
      for (size_t i = 1; i < args.size(); i++)
        if (args[i] < result)
          result = args[i];
      return result;
    }
    if (name == "max") {
      if (args.size() < 2)
        throw invalid_argument("max requires at least 2 arguments");
      BigInt result = args[0];
      for (size_t i = 1; i < args.size(); i++)
        if (args[i] > result)
          result = args[i];
      return result;
    }

    // Scientific functions (return integer approximation)
    if (name == "sin" || name == "cos" || name == "tan" || name == "asin" ||
        name == "acos" || name == "atan" || name == "sinh" || name == "cosh" ||
        name == "tanh" || name == "log" || name == "ln" || name == "log10" ||
        name == "exp") {
      if (args.size() != 1)
        throw invalid_argument(name + " requires 1 argument");
      double val = args[0].toDouble();
      double result;

      if (name == "sin")
        result = sin(calcState.useRadians ? val : toRadians(val));
      else if (name == "cos")
        result = cos(calcState.useRadians ? val : toRadians(val));
      else if (name == "tan")
        result = tan(calcState.useRadians ? val : toRadians(val));
      else if (name == "asin")
        result = calcState.useRadians ? asin(val) : toDegrees(asin(val));
      else if (name == "acos")
        result = calcState.useRadians ? acos(val) : toDegrees(acos(val));
      else if (name == "atan")
        result = calcState.useRadians ? atan(val) : toDegrees(atan(val));
      else if (name == "sinh")
        result = sinh(val);
      else if (name == "cosh")
        result = cosh(val);
      else if (name == "tanh")
        result = tanh(val);
      else if (name == "log" || name == "ln")
        result = log(val);
      else if (name == "log10")
        result = log10(val);
      else if (name == "exp")
        result = exp(val);
      else
        result = 0;

      return BigInt((long long)round(result));
    }

    throw invalid_argument("Unknown function: " + name);
  }

  BigInt parsePower() {
    BigInt left = parsePrimary();
    skip();
    if (peek() == '^') {
      get();
      BigInt right = parsePower(); // Right associative
      return power(left, right);
    }
    return left;
  }

  BigInt parseTerm() {
    BigInt left = parsePower();
    skip();
    while (peek() == '*' || peek() == '/' || peek() == '%') {
      char op = get();
      BigInt right = parsePower();
      if (op == '*')
        left *= right;
      else if (op == '/')
        left /= right;
      else
        left %= right;
      skip();
    }
    return left;
  }

  BigInt parseExpression() {
    skip();
    BigInt left = parseTerm();
    skip();
    while (peek() == '+' || peek() == '-') {
      char op = get();
      BigInt right = parseTerm();
      if (op == '+')
        left += right;
      else
        left -= right;
      skip();
    }
    return left;
  }

public:
  BigInt parse(const string &expression) {
    expr = expression;
    pos = 0;
    BigInt result = parseExpression();
    skip();
    if (pos != expr.size())
      throw invalid_argument("Unexpected characters at end");
    return result;
  }
};

// ═══════════════════════════════════════════════════════════════════════════
//  CALCULATOR COMMANDS
// ═══════════════════════════════════════════════════════════════════════════

void showHelp() {
  cout << R"(
  ╔═══════════════════════════════════════════════════════════════════════════╗
  ║                        ZILLION CALCULATOR HELP                            ║
  ╠═══════════════════════════════════════════════════════════════════════════╣
  ║  BASIC OPERATIONS                                                         ║
  ║    +, -, *, /         Addition, subtraction, multiplication, division    ║
  ║    %                  Modulo (remainder)                                 ║
  ║    ^                  Power (e.g., 2^10 = 1024)                          ║
  ║    ( )                Parentheses for grouping                           ║
  ╠═══════════════════════════════════════════════════════════════════════════╣
  ║  MATHEMATICAL FUNCTIONS                                                   ║
  ║    sqrt(x)            Square root                                        ║
  ║    cbrt(x), root(x,n) Cube root, nth root                               ║
  ║    abs(x)             Absolute value                                     ║
  ║    fact(n)            Factorial (n!)                                     ║
  ║    fib(n)             Nth Fibonacci number                               ║
  ║    pow(x,n)           Power function                                     ║
  ║    min(a,b,...), max(a,b,...)  Minimum/Maximum                          ║
  ╠═══════════════════════════════════════════════════════════════════════════╣
  ║  COMBINATORICS                                                            ║
  ║    nPr(n,r)           Permutation                                        ║
  ║    nCr(n,r)           Combination (binomial coefficient)                 ║
  ╠═══════════════════════════════════════════════════════════════════════════╣
  ║  NUMBER THEORY                                                            ║
  ║    gcd(a,b)           Greatest common divisor                            ║
  ║    lcm(a,b)           Least common multiple                              ║
  ║    prime(n)           Primality test (returns 1 if prime, 0 otherwise)   ║
  ╠═══════════════════════════════════════════════════════════════════════════╣
  ║  TRIGONOMETRIC (returns integer approximation)                            ║
  ║    sin(x), cos(x), tan(x)       Basic trig functions                     ║
  ║    asin(x), acos(x), atan(x)    Inverse trig functions                   ║
  ║    sinh(x), cosh(x), tanh(x)    Hyperbolic functions                     ║
  ╠═══════════════════════════════════════════════════════════════════════════╣
  ║  LOGARITHMIC                                                              ║
  ║    log(x), ln(x)      Natural logarithm                                  ║
  ║    log10(x)           Base-10 logarithm                                  ║
  ║    exp(x)             e^x                                                ║
  ╠═══════════════════════════════════════════════════════════════════════════╣
  ║  SPECIAL VARIABLES                                                        ║
  ║    ans                Previous result                                    ║
  ║    mem                Memory value                                       ║
  ║    pi                 Pi (integer: 3)                                    ║
  ╠═══════════════════════════════════════════════════════════════════════════╣
  ║  MEMORY COMMANDS                                                          ║
  ║    m+ or m+<expr>     Add to memory                                      ║
  ║    m- or m-<expr>     Subtract from memory                               ║
  ║    mr                 Memory recall (display memory)                     ║
  ║    mc                 Memory clear                                       ║
  ╠═══════════════════════════════════════════════════════════════════════════╣
  ║  OTHER COMMANDS                                                           ║
  ║    help               Show this help                                     ║
  ║    clear              Clear screen                                       ║
  ║    factors <n>        Show prime factorization                           ║
  ║    digits             Show digit count of last result                    ║
  ║    rad                Switch to radians mode                             ║
  ║    deg                Switch to degrees mode                             ║
  ║    exit, quit         Exit calculator                                    ║
  ╚═══════════════════════════════════════════════════════════════════════════╝
)";
}

void displayBanner() {
  cout << R"(
  ╔════════════════════════════════════════════════════════════╗
  ║                                                            ║
  ║    ███████╗██╗██╗     ██╗     ██╗ ██████╗ ███╗   ██╗       ║
  ║    ╚══███╔╝██║██║     ██║     ██║██╔═══██╗████╗  ██║       ║
  ║      ███╔╝ ██║██║     ██║     ██║██║   ██║██╔██╗ ██║       ║
  ║     ███╔╝  ██║██║     ██║     ██║██║   ██║██║╚██╗██║       ║
  ║    ███████╗██║███████╗███████╗██║╚██████╔╝██║ ╚████║       ║
  ║    ╚══════╝╚═╝╚══════╝╚══════╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝       ║
  ║                                                            ║
  ║       SCIENTIFIC CALCULATOR - Arbitrary Precision          ║
  ║                                                            ║
  ║  Type 'help' for commands, 'exit' to quit                  ║
  ╚════════════════════════════════════════════════════════════╝
)";
}

void processCommand(const string &input) {
  string cmd = trim(input);
  string lower = toLower(cmd);

  if (cmd.empty())
    return;

  // Commands
  if (lower == "help" || lower == "?") {
    showHelp();
    return;
  }

  if (lower == "exit" || lower == "quit" || lower == "q") {
    cout << "\n  Goodbye!\n\n";
    exit(0);
  }

  if (lower == "clear" || lower == "cls") {
    system(CLEAR_SCREEN);
    displayBanner();
    return;
  }

  if (lower == "mc") {
    calcState.memory = BigInt(0LL);
    cout << "  Memory cleared\n";
    return;
  }

  if (lower == "mr") {
    cout << "  Memory: " << calcState.memory << "\n";
    return;
  }

  if (lower == "rad") {
    calcState.useRadians = true;
    cout << "  Angle mode: Radians\n";
    return;
  }

  if (lower == "deg") {
    calcState.useRadians = false;
    cout << "  Angle mode: Degrees\n";
    return;
  }

  if (lower == "digits") {
    if (calcState.hasLastResult) {
      cout << "  Last result has " << calcState.lastResult.digitCount()
           << " digits\n";
    } else {
      cout << "  No previous result\n";
    }
    return;
  }

  if (lower.substr(0, 7) == "factors") {
    string numStr = trim(cmd.substr(7));
    if (numStr.empty() && calcState.hasLastResult) {
      numStr = calcState.lastResult.toString();
    }
    try {
      BigInt n(numStr);
      auto factors = primeFactorize(n);
      cout << "  Prime factorization of " << n << ":\n  ";
      for (size_t i = 0; i < factors.size(); i++) {
        if (i > 0)
          cout << " × ";
        cout << factors[i].first;
        if (factors[i].second > 1)
          cout << "^" << factors[i].second;
      }
      cout << "\n";
    } catch (...) {
      cout << "  Error: Invalid number\n";
    }
    return;
  }

  // Memory operations
  if (lower.substr(0, 2) == "m+") {
    string expr = trim(cmd.substr(2));
    BigInt val = calcState.hasLastResult ? calcState.lastResult : BigInt(0LL);
    if (!expr.empty()) {
      try {
        ExpressionParser parser;
        val = parser.parse(expr);
      } catch (...) {
      }
    }
    calcState.memory += val;
    cout << "  Memory: " << calcState.memory << "\n";
    return;
  }

  if (lower.substr(0, 2) == "m-") {
    string expr = trim(cmd.substr(2));
    BigInt val = calcState.hasLastResult ? calcState.lastResult : BigInt(0LL);
    if (!expr.empty()) {
      try {
        ExpressionParser parser;
        val = parser.parse(expr);
      } catch (...) {
      }
    }
    calcState.memory -= val;
    cout << "  Memory: " << calcState.memory << "\n";
    return;
  }

  // Evaluate expression
  try {
    ExpressionParser parser;
    BigInt result = parser.parse(cmd);
    calcState.lastResult = result;
    calcState.hasLastResult = true;

    cout << "  = " << result;
    if (result.digitCount() > 20) {
      cout << " (" << result.digitCount() << " digits)";
    }
    cout << "\n";
  } catch (const exception &e) {
    cout << "  Error: " << e.what() << "\n";
  }
}

// ═══════════════════════════════════════════════════════════════════════════
//  MAIN FUNCTION - REPL CALCULATOR
// ═══════════════════════════════════════════════════════════════════════════

int main() {
  system(CLEAR_SCREEN);
  displayBanner();

  string line;
  while (true) {
    cout << "\n  >> ";
    if (!getline(cin, line))
      break;
    processCommand(line);
  }

  cout << "\n  Goodbye!\n\n";
  return 0;
}
