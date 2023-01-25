#include "rational.h"

Rational::Rational() : numerator_(0), denominator_(1) {}

Rational::Rational(int num) : numerator_(num), denominator_(1) {}

Rational::Rational(const BigInteger& num) : numerator_(num), denominator_(1) {}

Rational::Rational(const BigInteger& numerator, const BigInteger& denominator)
    : numerator_(numerator), denominator_(denominator) {
    simplify();
}

Rational::Rational(const Rational& num)
    : numerator_(num.numerator_), denominator_(num.denominator_) {}

BigInteger Rational::getNumerator() const {
    return numerator_;
}

BigInteger Rational::getDenominator() const {
    return denominator_;
}

Rational& Rational::operator=(Rational num) {
    swap(num);
    return *this;
}

std::string Rational::toString() const {
    if (denominator_ == 1 || !numerator_) {
        return numerator_.toString();
    }
    return numerator_.toString() + "/" + denominator_.toString();
}

std::string Rational::asDecimal(size_t precision) const {
    auto [q, rem] = BigInteger::divide(numerator_, denominator_);
    std::string s;
    if (numerator_ < 0) {
        rem.changeSign();
        s += "-";
    }
    s += q.toString();
    if (precision == 0) {
        return s;
    }

    std::string after_dot;
    std::size_t end =
        (precision + BigInteger::base_len - 1) / BigInteger::base_len;
    for (size_t i = 0; i < end; ++i) {
        rem *= BigInteger::base;
        std::tie(q, rem) = BigInteger::divide(rem, denominator_);
        after_dot += q.toString(false);
    }
    after_dot = after_dot.substr(0, precision);
    return s + "." + after_dot;
}

Rational& Rational::operator+=(const Rational& num) {
    numerator_ = numerator_ * num.denominator_ + denominator_ * num.numerator_;
    denominator_ = denominator_ * num.denominator_;
    simplify();
    return *this;
}

Rational& Rational::operator-=(const Rational& num) {
    numerator_ = numerator_ * num.denominator_ - denominator_ * num.numerator_;
    denominator_ = denominator_ * num.denominator_;
    simplify();
    return *this;
}

Rational& Rational::operator*=(const Rational& num) {
    numerator_ *= num.numerator_;
    denominator_ *= num.denominator_;
    simplify();
    return *this;
}

Rational& Rational::operator/=(const Rational& num) {
    numerator_ *= num.denominator_;
    denominator_ *= num.numerator_;
    if (denominator_ < 0) {
        denominator_.changeSign();
        numerator_.changeSign();
    }
    simplify();
    return *this;
}

Rational Rational::operator-() const {
    return Rational(-numerator_, denominator_);
}

Rational operator+(Rational lhs, const Rational& rhs) {
    lhs += rhs;
    return lhs;
}

Rational operator-(Rational lhs, const Rational& rhs) {
    lhs -= rhs;
    return lhs;
}

Rational operator*(Rational lhs, const Rational& rhs) {
    lhs *= rhs;
    return lhs;
}

Rational operator/(Rational lhs, const Rational& rhs) {
    lhs /= rhs;
    return lhs;
}

bool operator==(const Rational& lhs, const Rational& rhs) {
    return lhs.getNumerator() == rhs.getNumerator() &&
           lhs.getDenominator() == rhs.getDenominator();
}

bool operator!=(const Rational& lhs, const Rational& rhs) {
    return !(lhs == rhs);
}

bool operator<(const Rational& lhs, const Rational& rhs) {
    return lhs.getNumerator() * rhs.getDenominator() <
           lhs.getDenominator() * rhs.getNumerator();
}

bool operator>(const Rational& lhs, const Rational& rhs) {
    return rhs < lhs;
}

bool operator<=(const Rational& lhs, const Rational& rhs) {
    return !(rhs < lhs);
}

bool operator>=(const Rational& lhs, const Rational& rhs) {
    return !(lhs < rhs);
}

std::istream& operator>>(std::istream& is, Rational& num) {
    BigInteger bi;
    is >> bi;
    num = Rational(bi);
    return is;
}

Rational::operator double() const {
    return std::stod(asDecimal(16));
}

void Rational::swap(Rational& num) {
    std::swap(numerator_, num.numerator_);
    std::swap(denominator_, num.denominator_);
}

void Rational::simplify() {
    if (!numerator_) {
        denominator_ = 1;
        return;
    }
    BigInteger d = BigInteger::gcd(numerator_ < 0 ? -numerator_ : numerator_,
                                   denominator_);
    numerator_ /= d;
    denominator_ /= d;
}
