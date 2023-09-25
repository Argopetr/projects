#include "CTreeSet.h"
#include "CBitSet.h"
#include "CListL2.h"
#include <ctime>
#include <cstdlib>

#define RND_DLT(x) (std::rand() % (2 * (x) + 1) - (x))
#define RND_DPS(a, b) ((a) + std::rand() % ((b) - (a) + 1))

#define SET(x) x

#define DELTA 100000

template <typename CSet>
void compair(const CSet&M1, const CSet&M2, const char* name1, const char* name2) {
    if ( M1 == M2 )
        std::cout << name1 << " == " << name2 << std::endl;
    else if ( M1 > M2 )
        std::cout << name1 << " > " << name2 << std::endl;
    else if ( M1 < M2 )
        std::cout << name1 << " < " << name2 << std::endl;
    else
        std::cout << name1 << " != " << name2 << std::endl;
}

template <typename CSet>
void ShowInfo(CSet& M, const char* name) {
    std::cout << name << " = " << M << std::endl;
    std::cout << "    AMOUNT OF ITEMS = " << M.abs() << std::endl;

    if ( M.isEmpty() )
        std::cout << "    NO ITEMS" << std::endl;
    else {
        int res = 0, tmp;

        for(tmp = M.begin(); M.end(); tmp = M.forward())
            res += tmp;

        std::cout << "    ITEMS SUM = " << res << std::endl;
        std::cout << "    MIN VAL = " << M.minVal() << std::endl;
        std::cout << "    MAX VAL = " << M.maxVal() << std::endl;
    }

    std::cout << std::endl;
}

void BitSetTest(void) {
    CBitSet M1, M2;
    int i, j;

    for(i = 0; i < RND_DPS(100, 1000); i++) {
        M1.makeEmpty();

        for(j = 0; j < RND_DPS(0, 1000); j++)
            M1.add(RND_DLT(DELTA));

        ShowInfo(M1, "M1");

        M2.makeEmpty();

        for(j = 0; j < RND_DPS(0, 1000); j++)
            M2.add(RND_DLT(DELTA));

        ShowInfo(M2, "M2");

        compair(M1, M2, "M1", "M2");
        std::cout << std::endl;

        std::cout << "M1 + M2 = " << M1 + M2 << std::endl;
        std::cout << "M1 - M2 = " << M1 - M2 << std::endl;
        std::cout << "M1 * M2 = " << M1 * M2 << std::endl;

        std::cout << "================================================================================" << std::endl << std::endl;
    }
}

void TreeSetTest(void) {
    CTreeSet<int> M1, M2;
    int i, j;

    for(i = 0; i < RND_DPS(100, 1000); i++) {
        M1.makeEmpty();

        for(j = 0; j < RND_DPS(0, 1000); j++)
            M1.add(RND_DLT(DELTA));

        ShowInfo(M1, "M1");

        M2.makeEmpty();

        for(j = 0; j < RND_DPS(0, 1000); j++)
            M2.add(RND_DLT(DELTA));

        ShowInfo(M2, "M2");

        compair(M1, M2, "M1", "M2");
        std::cout << std::endl;

        std::cout << "M1 + M2 = " << M1 + M2 << std::endl;
        std::cout << "M1 - M2 = " << M1 - M2 << std::endl;
        std::cout << "M1 * M2 = " << M1 * M2 << std::endl;

        std::cout << "================================================================================" << std::endl << std::endl;
    }
}

void ListL2Test(void) {
    CListL2<int> L1, L2, L3;
    int i, j;

    for(i = 0; i < RND_DPS(100, 1000); i++) {
        L1.clean();

        for(j = 0; j < RND_DPS(0, 1000); j++)
            L1.add(RND_DLT(DELTA));

        std::cout << "L1 = " << L1 << std::endl << std::endl;

        L2.clean();

        for(j = 0; j < RND_DPS(0, 1000); j++)
            L2.add(RND_DLT(DELTA));

        std::cout << "L2 = " << L2 << std::endl << std::endl;

        std::cout << "L1 join L2 = " << L1 + L2 << std::endl << std::endl;

	L3 = L1;

	for(i = 0; i < L2.amount(); i++)
		L3.del(L3.search(L2[i]));

	std::cout << "L1 without L2 elements = " << L3 << std::endl << std::endl;
        std::cout << "L1 sorted by increase = " << L1.sortInc() << std::endl << std::endl;
        std::cout << "L2 sorted by decrease = " << L2.sortDec() << std::endl << std::endl;

        std::cout << "================================================================================" << std::endl << std::endl;
    }
}

int main(void) {
    std::srand(std::time(0));
    BitSetTest();
    TreeSetTest();
    ListL2Test();

    return 0;
}
