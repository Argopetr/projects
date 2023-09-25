#include "Linal.h"
#include "SComplex.h"

int main() {
    SMatrix<SComplex> M;
    SVector<SComplex> v;

    std::cout << "INPUT SIZE AND ITEMS OF YOUR MATRIX M:" << std::endl;
    std::cin >> M;
    std::cout << std::endl;

    std::cout << "INPUT SIZE AND ITEMS OF YOUR VECTOR V:" << std::endl;
    std::cin >> v;
    std::cout << std::endl;

    std::cout << "PARTICULAR SOLUTION OF EQUATION MX = V:" << std::endl;
    std::cout << M.solve(v) << std::endl << std::endl;

    std::cout << "NULL SPACE" << std::endl;
    std::cout << M.nullSpace() << std::endl;

    return 0;
}
