#include <iostream>
#include "matrix.h"

int main() {
    const Matrix<5, 4, Residue<17>> bm = {{4, 0, 3, 2},
                                          {1, -7, 4, 5},
                                          {7, 1, 5, 3},
                                          {-5, -3, -3, -1},
                                          {1, -5, 2, 3}};
    std::cout << bm.rank() << "\n";
    std::cout << bm.det() << "\n";
    return 0;
}
