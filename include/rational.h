#include "biginteger.h"

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
