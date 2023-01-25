#include <algorithm>
#include <array>
#include "biginteger.h"
#include "rational.h"

namespace util {
template <size_t N, size_t D>
struct is_divisible {
    static const bool value = N % D == 0;
};

template <size_t N, size_t D>
const bool is_divisible_v = is_divisible<N, D>::value;

template <size_t N, size_t D>
struct is_square_less {
    static const bool value = D * D < N;
};

template <size_t N, size_t D>
const bool is_square_less_v = is_square_less<N, D>::value;

template <size_t N, size_t D, bool Condition>
struct is_not_prime {
    static const bool value =
        (is_divisible_v<N, D> ||
         is_not_prime<N, D + 1, is_square_less_v<N, D + 1>>::value);
};

template <size_t N, size_t D>
struct is_not_prime<N, D, false> {
    static const bool value = false;
};

template <size_t N, size_t D>
const bool is_not_prime_v = is_not_prime<N, D, is_square_less_v<N, D>>::value;

template <size_t N>
struct is_prime {
    static const bool value = !is_not_prime_v<N, 2>;
};

template <>
struct is_prime<2> {
    static const bool value = true;
};

template <size_t N>
const bool is_prime_v = is_prime<N>::value;
}  // namespace util

template <size_t N>
class Residue {
  public:
    Residue();
    explicit Residue(int x);
    explicit operator int();

    int getVal() const;
    Residue<N>& operator+=(const Residue<N>& another);
    Residue<N>& operator-=(const Residue<N>& another);
    Residue<N>& operator*=(const Residue<N>& another);
    Residue<N>& operator/=(const Residue<N>& another);

  private:
    int val;
};

template <size_t N>
Residue<N> operator+(Residue<N> lhs, const Residue<N>& rhs);

template <size_t N>
Residue<N> operator-(Residue<N> lhs, const Residue<N>& rhs);

template <size_t N>
Residue<N> operator*(Residue<N> lhs, const Residue<N>& rhs);

template <size_t N>
Residue<N> operator/(Residue<N> lhs, const Residue<N>& rhs);

template <size_t N>
bool operator==(const Residue<N>& lhs, const Residue<N>& rhs);

template <size_t N>
bool operator!=(const Residue<N>& lhs, const Residue<N>& rhs);

template <size_t N>
std::ostream& operator<<(std::ostream& os, const Residue<N>& elem);

template <size_t M, size_t N, typename Field = Rational>
class Matrix {
  public:
    Matrix();
    Matrix(std::initializer_list<std::initializer_list<int>> matrix);

    const std::array<Field, N>& getRow(size_t row) const;
    std::array<Field, M> getColumn(size_t col) const;

    Matrix<N, M, Field> transposed() const;
    Matrix<M, N, Field> getGaussMatrix() const;
    size_t rank() const;

    Field det() const;
    Field trace() const;
    void invert();
    Matrix<M, N, Field> inverted() const;

    Matrix<M, N, Field>& operator+=(const Matrix<M, N, Field>& another);
    Matrix<M, N, Field>& operator-=(const Matrix<M, N, Field>& another);
    Matrix<M, N, Field>& operator*=(const Field& val);
    Matrix<M, N, Field>& operator*=(const Matrix<M, N, Field>& another);
    std::array<Field, N>& operator[](size_t row);
    const std::array<Field, N>& operator[](size_t row) const;

  private:
    std::array<std::array<Field, N>, M> matrix_;
    int swaps_sign = 1;
};

template <size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator+(Matrix<M, N, Field> lhs,
                              const Matrix<M, N, Field>& rhs);

template <size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator-(Matrix<M, N, Field> lhs,
                              const Matrix<M, N, Field>& rhs);

template <size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator*(Matrix<M, N, Field> lhs, const Field& rhs);

template <size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator*(const Field& lhs, Matrix<M, N, Field> rhs);

template <size_t M, size_t N, size_t K, typename Field = Rational>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field>& lhs,
                              const Matrix<N, K, Field>& rhs);

template <size_t M, size_t N, typename Field = Rational>
bool operator==(const Matrix<M, N, Field>& lhs, const Matrix<M, N, Field>& rhs);

template <size_t M, size_t N, typename Field = Rational>
bool operator!=(const Matrix<M, N, Field>& lhs, const Matrix<M, N, Field>& rhs);

template <size_t M, typename Field = Rational>
using SquareMatrix = Matrix<M, M, Field>;

// Implementation

template <size_t N>
Residue<N>::Residue() : val(0) {}

template <size_t N>
Residue<N>::Residue(int x) : val(((x % static_cast<int>(N)) + N) % N) {}

template <size_t N>
Residue<N>::operator int() {
    return val;
}

template <size_t N>
int Residue<N>::getVal() const {
    return val;
}

template <size_t N>
Residue<N>& Residue<N>::operator+=(const Residue<N>& another) {
    val += another.val;
    val %= N;
    return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator-=(const Residue<N>& another) {
    val -= another.val - N;
    val %= N;
    return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator*=(const Residue<N>& another) {
    val *= another.val;
    val %= N;
    return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator/=(const Residue<N>& another) {
    static_assert(util::is_prime_v<N>, "N should be prime number in division!");
    Residue<N> inverse(1);
    for (size_t i = 0; i < N - 2; ++i) {
        inverse *= another;
    }
    val *= inverse.val;
    val %= N;
    return *this;
}

template <size_t N>
Residue<N> operator+(Residue<N> lhs, const Residue<N>& rhs) {
    return lhs += rhs;
}

template <size_t N>
Residue<N> operator-(Residue<N> lhs, const Residue<N>& rhs) {
    return lhs -= rhs;
}

template <size_t N>
Residue<N> operator*(Residue<N> lhs, const Residue<N>& rhs) {
    return lhs *= rhs;
}

template <size_t N>
Residue<N> operator/(Residue<N> lhs, const Residue<N>& rhs) {
    return lhs /= rhs;
}

template <size_t N>
bool operator==(const Residue<N>& lhs, const Residue<N>& rhs) {
    return lhs.getVal() == rhs.getVal();
}

template <size_t N>
bool operator!=(const Residue<N>& lhs, const Residue<N>& rhs) {
    return !(lhs == rhs);
}

template <size_t N>
std::ostream& operator<<(std::ostream& os, const Residue<N>& elem) {
    return os << elem.getVal();
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>::Matrix()
    : matrix_(std::array<std::array<Field, N>, M>()) {}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>::Matrix(
    std::initializer_list<std::initializer_list<int>> matrix) {
    size_t i = 0;
    for (const auto& r : matrix) {
        size_t j = 0;
        for (const int& v : r) {
            matrix_[i][j++] = Field(v);
        }
        ++i;
    }
}

template <size_t M, size_t N, typename Field>
const std::array<Field, N>& Matrix<M, N, Field>::getRow(size_t row) const {
    return matrix_[row];
}

template <size_t M, size_t N, typename Field>
std::array<Field, M> Matrix<M, N, Field>::getColumn(size_t col) const {
    std::array<Field, M> column;
    for (size_t i = 0; i < M; ++i) {
        column[i] = matrix_[i][col];
    }
    return column;
}

template <size_t M, size_t N, typename Field>
Matrix<N, M, Field> Matrix<M, N, Field>::transposed() const {
    Matrix<N, M, Field> transposed_matrix;
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            transposed_matrix[j][i] = matrix_[i][j];
        }
    }
    return transposed_matrix;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::getGaussMatrix() const {
    Matrix<M, N, Field> gauss_matrix = *this;
    size_t cur_row = 0;
    for (size_t col = 0; col < N; ++col) {
        bool swapped = false;
        for (size_t row = cur_row; row < M; ++row) {
            if (gauss_matrix[row][col] != Field(0)) {
                std::swap(gauss_matrix[cur_row], gauss_matrix[row]);
                swapped = true;
                if (cur_row != row) {
                    gauss_matrix.swaps_sign *= -1;
                }
                break;
            }
        }
        if (!swapped) {
            continue;
        }
        Field leader = gauss_matrix[cur_row][col];
        for (size_t row = cur_row + 1; row < M; ++row) {
            Field coefficient = gauss_matrix[row][col] / leader;
            if (coefficient == Field(0)) {
                continue;
            }
            for (size_t i = 0; i < N; ++i) {
                gauss_matrix[row][i] -= coefficient * gauss_matrix[cur_row][i];
            }
        }
        ++cur_row;
    }
    return gauss_matrix;
}

template <size_t M, size_t N, typename Field>
size_t Matrix<M, N, Field>::rank() const {
    Matrix<M, N, Field> gauss_matrix = getGaussMatrix();
    size_t rank = M;
    for (size_t row = M; row-- > 0;) {
        if (std::any_of(gauss_matrix[row].begin(), gauss_matrix[row].end(),
                        [](Field elem) {
                            return elem != Field(0);
                        })) {
            return rank;
        }
        --rank;
    }
    return rank;
}

template <size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::det() const {
    Matrix<M, N, Field> gauss_matrix = getGaussMatrix();
    Field d = Field(gauss_matrix.swaps_sign);
    for (size_t i = 0; i < M; ++i) {
        d *= gauss_matrix[i][i];
    }
    return d;
}

template <size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::trace() const {
    static_assert(M == N, "Matrix must be square!");
    Field tr = Field(0);
    for (size_t i = 0; i < M; ++i) {
        tr += matrix_[i][i];
    }
    return tr;
}

template <size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::invert() {
    static_assert(M == N, "Matrix must be square!");
    *this = inverted();
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::inverted() const {
    static_assert(M == N, "Matrix must be square!");
    Matrix<M, M, Field> inv;
    for (size_t i = 0; i < M; ++i) {
        inv[i][i] = Field(1);
    }

    Matrix<M, M, Field> gauss_matrix = *this;
    size_t cur_row = 0;
    for (size_t col = 0; col < M; ++col) {
        bool swapped = false;
        for (size_t row = cur_row; row < M; ++row) {
            if (gauss_matrix[row][col] != Field(0)) {
                std::swap(gauss_matrix[cur_row], gauss_matrix[row]);
                std::swap(inv[cur_row], inv[row]);
                swapped = true;
                break;
            }
        }
        if (!swapped) {
            continue;
        }
        Field leader = gauss_matrix[cur_row][col];
        for (size_t row = cur_row + 1; row < M; ++row) {
            Field coefficient = gauss_matrix[row][col] / leader;
            if (coefficient == Field(0)) {
                continue;
            }
            for (size_t i = 0; i < M; ++i) {
                gauss_matrix[row][i] -= coefficient * gauss_matrix[cur_row][i];
                inv[row][i] -= coefficient * inv[cur_row][i];
            }
        }
        ++cur_row;
    }

    for (size_t i = M; i-- > 0;) {
        Field leader = gauss_matrix[i][i];
        for (size_t row = 0; row < i; ++row) {
            Field coefficient = gauss_matrix[row][i] / leader;
            if (coefficient == Field(0)) {
                continue;
            }
            for (size_t j = 0; j < M; ++j) {
                gauss_matrix[row][j] -= coefficient * gauss_matrix[i][j];
                inv[row][j] -= coefficient * inv[i][j];
            }
        }
        for (size_t j = 0; j < M; ++j) {
            inv[i][j] /= gauss_matrix[i][i];
        }
    }

    return inv;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator+=(
    const Matrix<M, N, Field>& another) {
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            matrix_[i][j] += another[i][j];
        }
    }
    return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator-=(
    const Matrix<M, N, Field>& another) {
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            matrix_[i][j] -= another[i][j];
        }
    }
    return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator*=(const Field& val) {
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            matrix_[i][j] *= val;
        }
    }
    return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator*=(
    const Matrix<M, N, Field>& another) {
    static_assert(M == N, "Matrix must be square!");
    return *this = *this * another;
}

template <size_t M, size_t N, typename Field>
std::array<Field, N>& Matrix<M, N, Field>::operator[](size_t row) {
    return matrix_[row];
}

template <size_t M, size_t N, typename Field>
const std::array<Field, N>& Matrix<M, N, Field>::operator[](size_t row) const {
    return matrix_[row];
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator+(Matrix<M, N, Field> lhs,
                              const Matrix<M, N, Field>& rhs) {
    return lhs += rhs;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator-(Matrix<M, N, Field> lhs,
                              const Matrix<M, N, Field>& rhs) {
    return lhs -= rhs;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(Matrix<M, N, Field> lhs, const Field& rhs) {
    return lhs *= rhs;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Field& lhs, Matrix<M, N, Field> rhs) {
    return rhs *= lhs;
}

template <size_t M, size_t N, size_t K, typename Field>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field>& lhs,
                              const Matrix<N, K, Field>& rhs) {
    Matrix<M, K, Field> res;
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < K; ++j) {
            for (size_t k = 0; k < N; ++k) {
                res[i][j] += lhs[i][k] * rhs[k][j];
            }
        }
    }
    return res;
}

template <size_t M, size_t N, typename Field>
bool operator==(const Matrix<M, N, Field>& lhs,
                const Matrix<M, N, Field>& rhs) {
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (lhs[i][j] != rhs[i][j]) {
                return false;
            }
        }
    }
    return true;
}

template <size_t M, size_t N, typename Field>
bool operator!=(const Matrix<M, N, Field>& lhs,
                const Matrix<M, N, Field>& rhs) {
    return !(lhs == rhs);
}
