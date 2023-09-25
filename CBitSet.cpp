#include "CBitSet.h"

#define MIN(a, b) (((a) <= (b))? (a) : (b))
#define MAX(a, b) (((a) >= (b))? (a) : (b))

#define INT_BIT sizeof(int)

#define GET(arr, ind) ((arr)[(ind) / INT_BIT] & (1 << (ind) % INT_BIT))
#define TURNON(arr, ind) (arr)[(ind) / INT_BIT] |= 1 << (ind) % INT_BIT
#define TURNOFF(arr, ind) (arr)[(ind) / INT_BIT] &= ~(1 << (ind) % INT_BIT)
#define INRANGE(set, x) (((x) >= (set).StartValue && (x) < (set).StartValue + (set).Lenght * INT_BIT)? 1 : 0)

CBitSet::CBitSet() {
    StartValue = 0;
    Lenght = 0;
    Crr = 0;

    LBorder = 0;
    RBorder = -1;

    Space = NULL;
}

CBitSet::CBitSet(const CBitSet& M) {
    StartValue = M.StartValue;
    Lenght = M.Lenght;
    Crr = M.Crr;

    LBorder = M.LBorder;
    RBorder = M.RBorder;

    Space = new unsigned int[M.Lenght];

    for(int i = 0; i < M.Lenght; i++)
        Space[i] = M.Space[i];
}

CBitSet::CBitSet(int a, int b) {
    StartValue = a;
    Crr = 0;

    LBorder = a;
    RBorder = b;

    if ( a <= b ) {
        Lenght = 1;
        StartValue = a;

        for(; Lenght * INT_BIT < b - a; Lenght <<= 1);

        Space = new unsigned int[Lenght];
        makeEmpty();
    } else {
        Lenght = 0;
        Space = NULL;
    }
}

CBitSet& CBitSet::makeEmpty(void) {
    for(int i = 0; i < Lenght; i++)
        Space[i] = 0;

    return *this;
}

CBitSet& CBitSet::makeEmpty(int a, int b) {
    delete[] Space;

    StartValue = a;
    Crr = 0;

    LBorder = a;
    RBorder = b;

    if ( a <= b ) {
        Lenght = 1;
        StartValue = a;

        for(; Lenght * INT_BIT < b - a; Lenght <<= 1);

        Space = new unsigned int[Lenght];
        makeEmpty();
    } else {
        Lenght = 0;
        Space = NULL;
    }

    return *this;
}

int CBitSet::isEmpty(void) const {
    int res = 1;

    for(int i = 0; res == 1 && i < Lenght; i++)
        if ( Space[i] )
            res = 0;

    return res;
}

CBitSet& CBitSet::cut(void) {
	int a = 0, b = Lenght * INT_BIT - 1;

    for(; a < Lenght * INT_BIT && !GET(Space, a); a++);
    for(; b >= 0 && !GET(Space, b); b--);

    LBorder = a;
    RBorder = b;

    if ( a <= b ) {
        unsigned int *newSpace;
        int len = 1;

        for(; len * INT_BIT < b - a + 1; len++);

        newSpace = new unsigned int[len];

        for(int i = 0; i < b - a + 1; i++)
            if ( GET(Space, i + a) )
                TURNON(newSpace, i);
            else
                TURNOFF(newSpace, i);

        for(int i = b - a + 1; i < len * INT_BIT; i++)
            TURNOFF(newSpace, i);

        StartValue += a;
        Lenght = len;

        delete[] Space;
        Space = newSpace;
    } else {
        delete[] Space;

        Lenght = 0;
        Space = NULL;
    }

    return *this;
}

CBitSet& CBitSet::add(int val) {
    if ( !Lenght ) {
        StartValue = val;
        Lenght = 1;

        LBorder = val;
        RBorder = val;

        Space = new unsigned int[1];
        Space[0] = 1;
    } else if ( val < StartValue ) {
        int len = Lenght * 2;
        unsigned int *newSpace;

        LBorder = val;

        for(; (len - Lenght) * INT_BIT < StartValue - val; len *= 2);

        newSpace = new unsigned int[len];
        newSpace[0] = 1;

        for(int i = 1; i < len; i++)
            newSpace[i] = 0;

        for(int i = 0; i < Lenght * INT_BIT; i++)
            if ( GET(Space, i) )
                TURNON(newSpace, i + StartValue - val);

        delete[] Space;

        StartValue = val;
        Lenght = len;
        Space = newSpace;
    } else if ( val >= StartValue + Lenght * INT_BIT ) {
        int len = Lenght * 2;
        unsigned int *newSpace;

        RBorder = val;

        for(; (len - Lenght) * INT_BIT < val - StartValue - Lenght * INT_BIT; len *= 2);

        newSpace = new unsigned int[len];

        for(int i = 0; i < len; i++)
            newSpace[i] = 0;

        TURNON(newSpace, val - StartValue);

        for(int i = 0; i < Lenght * INT_BIT; i++)
            if ( GET(Space, i) )
                TURNON(newSpace, i);

        delete[] Space;

        Lenght = len;
        Space = newSpace;
    } else
        TURNON(Space, val - StartValue);

    return *this;
}

CBitSet& CBitSet::del(int val) {
    if ( INRANGE(*this, val) )
        TURNOFF(Space, val - StartValue);

    return *this;
}

int CBitSet::contains(int val) const {
    return (INRANGE(*this, val) && GET(Space, val - StartValue))? 1 : 0;
}

CBitSet CBitSet::Union(const CBitSet& M) const {
    CBitSet res(MIN(this -> LBorder, M.LBorder), MAX(this -> RBorder, M.RBorder));

    for(int i = 0; i < Lenght * INT_BIT; i++)
        if ( GET(Space, i) )
            TURNON(res.Space, i + this -> LBorder - res.LBorder);

    for(int j = 0; j < M.Lenght * INT_BIT; j++)
        if ( GET(M.Space, j) )
            TURNON(res.Space, j + M.leftBorder() - res.leftBorder());

    res.cut();

    return res;
}

CBitSet CBitSet::Intersection(const CBitSet& M) const {
    CBitSet res(MAX(this -> LBorder, M.LBorder), MIN(this -> RBorder, M.RBorder));

    for(int i = 0; i < res.Lenght * INT_BIT; i++)
        if ( this -> contains(i + res.StartValue) && M.contains(i + res.StartValue) )
            TURNON(res.Space, i);

    res.cut();

    return res;
}

CBitSet CBitSet::Difference(const CBitSet& M) const {
    CBitSet res(*this);

    for(int i = 0; i < Lenght * INT_BIT; i++)
        if ( M.contains(i + StartValue) == 1 )
            TURNOFF(res.Space, i);

    res.cut();

    return res;
}

int CBitSet::minVal(void) const {
    int res = 0;

    for(; res < Lenght * INT_BIT && !GET(Space, res); res++);

    res += StartValue;

    return res;

}

int CBitSet::maxVal(void) const {
    int res = Lenght * INT_BIT - 1;

    for(; res >= 0 && !GET(Space, res); res--);

    res += StartValue;

    return res;
}

int CBitSet::abs(void) const {
    int res = 0, i;

    for(i = 0; i < Lenght * INT_BIT; i++)
        if ( GET(Space, i) )
            res++;

    return res;
}

int CBitSet::begin(void) {
    Crr = 0;

    for(; Crr < Lenght * INT_BIT && !GET(Space, Crr); Crr++);

    return Crr + StartValue;
}

int CBitSet::forward(void) {
    for( Crr++; Crr < Lenght * INT_BIT && !GET(Space, Crr); Crr++);

    return Crr + StartValue;
}

int CBitSet::end(void) const {
    return (Crr >= Lenght * INT_BIT)? 0 : 1;
}

CBitSet& CBitSet::operator=(const CBitSet& M) {
    if ( &M != this ) {
        StartValue = M.StartValue;
        Lenght = M.Lenght;
        Crr = M.Crr;

        LBorder = M.LBorder;
        RBorder = M.RBorder;

        delete[] Space;
        Space = new unsigned int[M.Lenght];

        for(int i = 0; i < M.Lenght; i++)
            Space[i] = M.Space[i];
    }

    return *this;
}

int CBitSet::operator<=(const CBitSet& M) const {
    int res = 1;

    for(int i = 0; res && i < Lenght; i++)
        if ( GET(Space, i) && !M.contains(i + StartValue - M.StartValue) )
            res = 0;

    return res;
}

CSetBuffer<CBitSet, int> CBitSet::operator[](const int& x) {
    CSetBuffer<CBitSet, int> res(this, x);

    return res;
}

std::istream& operator>>(std::istream& input, CBitSet& M) {
    int num, tmp;

    input >> num;

    for(; num > 0; num--) {
        input >> tmp;
        M.add(tmp);
    }

    return input;
}

std::ostream& operator<<(std::ostream& output, const CBitSet& M) {
    if ( M.isEmpty() )
        output << "Empty set";
    else {
        int flag = 0;

        output << '{';

        for(int i = 0; i < M.Lenght * INT_BIT; i++)
            if ( GET(M.Space, i) ) {
                if ( flag )
                    output << ", ";
                else
                    flag = 1;

                output << i + M.StartValue;
            }

        output << '}';
    }

    return output;
}
