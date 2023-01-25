#pragma once

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

class BigInteger {
  public:
    static const long long base_len = 9;
    static const long long base = 1e9;

    static std::pair<BigInteger, BigInteger> divide(const BigInteger& n,
                                                    const BigInteger& d);
    static BigInteger gcd(const BigInteger& a, const BigInteger& b);

    BigInteger();
    BigInteger(int x);
    BigInteger(std::string str);
    BigInteger(const BigInteger& num);

    std::string toString(bool no_leading_zeros = true) const;
    int compare(const BigInteger& num) const;
    void changeSign();
    void shift();

    BigInteger& operator=(BigInteger num);
    BigInteger& operator+=(const BigInteger& num);
    BigInteger& operator-=(const BigInteger& num);
    BigInteger& operator*=(const BigInteger& num);
    BigInteger& operator/=(const BigInteger& num);
    BigInteger& operator%=(const BigInteger& num);
    BigInteger operator-() const;
    BigInteger& operator++();
    BigInteger operator++(int);
    BigInteger& operator--();
    BigInteger operator--(int);

    explicit operator bool() const;

  private:
    bool is_negative_ = false;
    std::vector<long long> digits_;

    static size_t getNumberLength(long long x);
    static size_t getNumberLength(std::string str);
    static long long findBaseDivider(const BigInteger& n, const BigInteger& d);

    void swap(BigInteger& num);
    void removeLeadingZeros();
};

BigInteger operator+(BigInteger lhs, const BigInteger& rhs);
BigInteger operator-(BigInteger lhs, const BigInteger& rhs);
BigInteger operator*(BigInteger lhs, const BigInteger& rhs);
BigInteger operator/(BigInteger lhs, const BigInteger& rhs);
BigInteger operator%(BigInteger lhs, const BigInteger& rhs);
bool operator==(const BigInteger& lhs, const BigInteger& rhs);
bool operator!=(const BigInteger& lhs, const BigInteger& rhs);
bool operator<(const BigInteger& lhs, const BigInteger& rhs);
bool operator>(const BigInteger& lhs, const BigInteger& rhs);
bool operator<=(const BigInteger& lhs, const BigInteger& rhs);
bool operator>=(const BigInteger& lhs, const BigInteger& rhs);
BigInteger operator""_bi(unsigned long long x);
BigInteger operator""_bi(const char* str, size_t /*unused*/);

std::ostream& operator<<(std::ostream& os, const BigInteger& num);
std::istream& operator>>(std::istream& is, BigInteger& num);

class Rational {
  public:
    Rational();
    Rational(int num);
    Rational(const BigInteger& num);
    Rational(const BigInteger& numerator, const BigInteger& denominator);
    Rational(const Rational& num);

    BigInteger getNumerator() const;
    BigInteger getDenominator() const;

    std::string toString() const;
    std::string asDecimal(size_t precision = 0) const;

    Rational& operator=(Rational num);
    Rational& operator+=(const Rational& num);
    Rational& operator-=(const Rational& num);
    Rational& operator*=(const Rational& num);
    Rational& operator/=(const Rational& num);
    Rational operator-() const;

    explicit operator double() const;

  private:
    BigInteger numerator_;
    BigInteger denominator_;

    void swap(Rational& num);
    void simplify();
};

Rational operator+(Rational lhs, const Rational& rhs);
Rational operator-(Rational lhs, const Rational& rhs);
Rational operator*(Rational lhs, const Rational& rhs);
Rational operator/(Rational lhs, const Rational& rhs);
bool operator==(const Rational& lhs, const Rational& rhs);
bool operator!=(const Rational& lhs, const Rational& rhs);
bool operator<(const Rational& lhs, const Rational& rhs);
bool operator>(const Rational& lhs, const Rational& rhs);
bool operator<=(const Rational& lhs, const Rational& rhs);
bool operator>=(const Rational& lhs, const Rational& rhs);
std::istream& operator>>(std::istream& is, Rational& num);
