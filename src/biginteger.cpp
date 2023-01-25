#include "biginteger.h"

size_t BigInteger::getNumberLength(long long x) {
    if (x == 0) {
        return 1;
    }
    return std::ceil(std::log(std::abs(x) + 1) / std::log(base));
}

size_t BigInteger::getNumberLength(std::string str) {
    size_t sz = str.size();
    if (str[0] == '-') {
        --sz;
    }
    return (sz + base_len - 1) / base_len;
}

BigInteger::BigInteger() : is_negative_(false) {}

BigInteger::BigInteger(int x)
    : is_negative_(x < 0), digits_(getNumberLength(x), 0) {
    x = std::abs(x);
    for (size_t i = 0; i < digits_.size(); ++i) {
        digits_[i] = x % base;
        x /= base;
    }
}

BigInteger::BigInteger(std::string str)
    : is_negative_(str[0] == '-'), digits_(getNumberLength(str), 0) {
    size_t k = 0;
    int start = str.size();
    int end = 0;
    if (is_negative_) {
        ++end;
    }
    for (int i = start; i > end; i -= base_len) {
        if (i >= base_len) {
            digits_[k++] = stoi(str.substr(i - base_len, base_len));
        } else {
            digits_[k++] = stoi(str.substr(end, i - end));
        }
    }
    removeLeadingZeros();
}

BigInteger::BigInteger(const BigInteger& num)
    : is_negative_(num.is_negative_), digits_(num.digits_) {}

std::string BigInteger::toString(bool no_leading_zeros) const {
    std::stringstream ss;
    if (is_negative_) {
        ss << '-';
    }
    auto iter = digits_.rbegin();
    ss << std::setfill('0');
    if (!no_leading_zeros) {
        ss << std::setw(base_len);
    }
    ss << *iter++;
    for (; iter != digits_.rend(); ++iter) {
        ss << std::setw(base_len) << *iter;
    }
    return ss.str();
}

void BigInteger::changeSign() {
    if (*this != 0) {
        is_negative_ ^= 1;
    }
}

BigInteger BigInteger::gcd(const BigInteger& a, const BigInteger& b) {
    if (b == 0) {
        return a;
    }
    return gcd(b, a % b);
}

void BigInteger::shift() {
    if (*this != 0) {
        digits_.insert(digits_.begin(), 0);
    }
}

BigInteger& BigInteger::operator=(BigInteger num) {
    swap(num);
    return *this;
}

BigInteger& BigInteger::operator+=(const BigInteger& num) {
    if (is_negative_ != num.is_negative_) {
        if (is_negative_) {
            return *this = num - (-*this);
        }
        return *this -= -num;
    }

    long long carry = 0;
    size_t n = digits_.size();
    size_t m = num.digits_.size();
    size_t max_size = std::max(n, m);
    for (size_t i = 0; i < max_size || carry != 0; ++i) {
        long long x = carry;
        if (i < n) {
            x += digits_[i];
        }
        if (i < m) {
            x += num.digits_[i];
        }
        carry = static_cast<long long>(x >= base);
        if (carry != 0) {
            x -= base;
        }

        if (i >= n) {
            digits_.push_back(0);
        }
        digits_[i] = x;
    }

    return *this;
}

BigInteger& BigInteger::operator-=(const BigInteger& num) {
    if (num.is_negative_) {
        return *this += -num;
    }
    if (*this == num) {
        return *this = 0;
    }
    if (*this < num) {
        *this = num - *this;
        is_negative_ = true;
        return *this;
    }

    size_t n = digits_.size();
    size_t m = num.digits_.size();
    long long carry = 0;
    for (size_t i = 0; i < n; ++i) {
        long long x = -carry;
        if (i < n) {
            x += digits_[i];
        }
        if (i < m) {
            x -= num.digits_[i];
        }
        carry = static_cast<long long>(x < 0);
        if (carry != 0) {
            x += base;
        }
        digits_[i] = x;
    }

    removeLeadingZeros();
    return *this;
}

BigInteger& BigInteger::operator*=(const BigInteger& num) {
    size_t n = digits_.size();
    size_t m = num.digits_.size();
    if (n < m) {
        return *this = num * (*this);
    }

    std::vector<long long> new_digits(n + m);
    for (size_t i = 0; i < n; ++i) {
        long long carry = 0;
        for (size_t j = 0; j < m || carry != 0; ++j) {
            long long x = new_digits[i + j] + carry;
            if (j < m) {
                x += digits_[i] * num.digits_[j];
            }
            new_digits[i + j] = x % base;
            carry = x / base;
        }
    }

    digits_ = std::move(new_digits);
    removeLeadingZeros();
    is_negative_ = digits_.back() != 0 && (is_negative_ != num.is_negative_);
    return *this;
}

std::pair<BigInteger, BigInteger> BigInteger::divide(const BigInteger& n,
                                                     const BigInteger& d) {
    if (d.is_negative_) {
        auto p = divide(n, -d);
        return {-p.first, p.second};
    }
    if (n.is_negative_) {
        auto p = divide(-n, d);
        return {-p.first, -p.second};
    }

    if (n < d) {
        return {0, n};
    }
    BigInteger q = 0;
    BigInteger r = 0;
    for (size_t i = n.digits_.size(); i-- > 0;) {
        r.shift();
        r += n.digits_[i];
        if (r < d) {
            q.shift();
            continue;
        }
        long long j = findBaseDivider(r, d);
        r -= d * j;
        q.shift();
        q += j;
    }

    return {q, r};
}

BigInteger& BigInteger::operator/=(const BigInteger& num) {
    return *this = divide(*this, num).first;
}

BigInteger& BigInteger::operator%=(const BigInteger& num) {
    return *this = divide(*this, num).second;
}

BigInteger BigInteger::operator-() const {
    BigInteger num(*this);
    if (num != 0) {
        num.is_negative_ ^= 1;
    }
    return num;
}

BigInteger& BigInteger::operator++() {
    return *this += 1;
}

BigInteger BigInteger::operator++(int) {
    BigInteger copy = *this;
    ++*this;
    return copy;
}

BigInteger& BigInteger::operator--() {
    return *this -= 1;
}

BigInteger BigInteger::operator--(int) {
    BigInteger copy = *this;
    --*this;
    return copy;
}

BigInteger operator+(BigInteger lhs, const BigInteger& rhs) {
    return lhs += rhs;
}

BigInteger operator-(BigInteger lhs, const BigInteger& rhs) {
    return lhs -= rhs;
}

BigInteger operator*(BigInteger lhs, const BigInteger& rhs) {
    return lhs *= rhs;
}

BigInteger operator/(BigInteger lhs, const BigInteger& rhs) {
    return lhs /= rhs;
}

BigInteger operator%(BigInteger lhs, const BigInteger& rhs) {
    return lhs %= rhs;
}

// -1 - less.
// 0 - equal.
// 1 - greater.
int BigInteger::compare(const BigInteger& num) const {
    if (is_negative_ != num.is_negative_) {
        return is_negative_ ? -1 : 1;
    }
    size_t n = digits_.size();
    size_t m = num.digits_.size();
    if (n != m) {
        return (is_negative_ ^ (n < m)) != 0 ? -1 : 1;
    }
    for (size_t i = n; i-- > 0;) {
        if (digits_[i] != num.digits_[i]) {
            return (is_negative_ ^ (digits_[i] < num.digits_[i])) != 0 ? -1 : 1;
        }
    }
    return 0;
}

bool operator==(const BigInteger& lhs, const BigInteger& rhs) {
    return lhs.compare(rhs) == 0;
}

bool operator!=(const BigInteger& lhs, const BigInteger& rhs) {
    return !(lhs == rhs);
}

bool operator<(const BigInteger& lhs, const BigInteger& rhs) {
    return lhs.compare(rhs) == -1;
}

bool operator>(const BigInteger& lhs, const BigInteger& rhs) {
    return rhs < lhs;
}

bool operator<=(const BigInteger& lhs, const BigInteger& rhs) {
    return !(rhs < lhs);
}

bool operator>=(const BigInteger& lhs, const BigInteger& rhs) {
    return !(lhs < rhs);
}

BigInteger::operator bool() const {
    return digits_.back() != 0;
}

BigInteger operator""_bi(unsigned long long x) {
    return BigInteger(std::to_string(x));
}

BigInteger operator""_bi(const char* str, size_t /*unused*/) {
    return BigInteger(std::string(str));
}

std::ostream& operator<<(std::ostream& os, const BigInteger& num) {
    return os << num.toString();
}

std::istream& operator>>(std::istream& is, BigInteger& num) {
    std::string s;
    is >> s;
    num = s;
    return is;
}

void BigInteger::swap(BigInteger& num) {
    std::swap(is_negative_, num.is_negative_);
    std::swap(digits_, num.digits_);
}

void BigInteger::removeLeadingZeros() {
    while (digits_.size() > 1 && digits_.back() == 0) {
        digits_.pop_back();
    }
}

long long BigInteger::findBaseDivider(const BigInteger& n,
                                      const BigInteger& d) {
    long long lo = 1, hi = base;
    while (hi - lo > 1) {
        long long mid = lo + (hi - lo) / 2;
        if (mid * d <= n) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}
