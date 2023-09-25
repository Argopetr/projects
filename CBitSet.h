#ifndef CBITSET_H
#define CBITSET_H

#include <iostream>
#include "CSetBuffer.h"

struct CBitSet {
private:
    int StartValue, Lenght, Crr;
    int LBorder, RBorder;
    unsigned int *Space;

public:
    CBitSet();
    CBitSet(const CBitSet& M);
    CBitSet(int a, int b);
    ~CBitSet() {delete[] Space;}

    CBitSet& makeEmpty(void);
    CBitSet& makeEmpty(int a, int b);
    int isEmpty(void) const;

    CBitSet& makeInvalid(void) {return this -> makeEmpty(0, -1);}
    int isValid(void) const {return (Lenght)? 1 : 0;}

    CBitSet& cut(void);
    CBitSet& add(int val);
    CBitSet& del(int val);
    int contains(int val) const;

    CBitSet Intersection(const CBitSet& M) const;
    CBitSet Union(const CBitSet& M) const;
    CBitSet Difference(const CBitSet& M) const;

    int leftBorder(void) const {return LBorder;}
    int rightBorder(void) const {return RBorder;}

    int minVal(void) const;
    int maxVal(void) const;
    int abs(void) const;

    int begin(void); 
    int forward(void);
    int end(void) const;

    int operator*(int x) const {return this -> contains(x);}

    CBitSet& operator=(const CBitSet& M);

    CBitSet& operator*=(const CBitSet& M) {return *this = this -> Intersection(M);}
    CBitSet& operator+=(const CBitSet& M) {return *this = this -> Union(M);}
    CBitSet& operator-=(const CBitSet& M) {return *this = this -> Difference(M);}

    CBitSet operator*(const CBitSet& M) {return this -> Intersection(M);}
    CBitSet operator+(const CBitSet& M) {return this -> Union(M);}
    CBitSet operator-(const CBitSet& M) {return this -> Difference(M);}

    CBitSet& operator+=(int x) {return this -> add(x);}
    CBitSet& operator-=(int x) {return this -> del(x);}

    CBitSet operator+(int x) {CBitSet res(*this); return res.add(x);}
    CBitSet operator-(int x) {CBitSet res(*this); return res.del(x);}

    int operator<=(const CBitSet& M) const;
    int operator>=(const CBitSet& M) const {return M <= *this;}
    int operator==(const CBitSet& M) const {return (M <= *this) * (M >= *this);}
    int operator!=(const CBitSet& M) const {return 1 - (M == *this);}
    int operator<(const CBitSet& M) const {return (*this <= M) * ((this -> abs() < M.abs())? 1 : 0);}
    int operator>(const CBitSet& M) const {return M < *this;}

    CSetBuffer<CBitSet, int> operator[](const int& x);

    friend std::istream& operator>>(std::istream& input, CBitSet& M);
    friend std::ostream& operator<<(std::ostream& output, const CBitSet& M);
};

#endif
