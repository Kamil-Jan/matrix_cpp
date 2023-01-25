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
