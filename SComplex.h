#ifndef COMPLEX_H
#define COMPLEX_H

#include <iostream>
#include "Neutral.h"

struct SComplex {
private:
    double real, image;

public:
    SComplex() {real = image = 0;}
    SComplex(double x);
    SComplex(double a, double b);

    double& re(void) {return real;}
    double& im(void) {return image;}

    double abs(void) const {return std::sqrt(real * real + image * image);}
    double arg(void) const;

    SComplex& operator=(const double& x);
    SComplex& operator+=(const double& x);
    SComplex& operator-=(const double& x);
    SComplex& operator*=(const double& x);
    SComplex& operator/=(const double& x);

    SComplex operator+(const double& x) const;
    SComplex operator-(const double& x) const;
    SComplex operator*(const double& x) const;
    SComplex operator/(const double& x) const;

    SComplex& operator+=(const SComplex& x);
    SComplex& operator-=(const SComplex& x);
    SComplex& operator*=(const SComplex& x);
    SComplex& operator/=(const SComplex& x);

    SComplex operator+(const SComplex& x) const;
    SComplex operator-(const SComplex& x) const;
    SComplex operator*(const SComplex& x) const;
    SComplex operator/(const SComplex& x) const;

    SComplex& operator~(void);
    SComplex& operator=(double& x);

    char operator==(const double& x) const;
    char operator!=(const double& x) const {return 1 - ((*this) == x);}

    char operator==(const SComplex& x) const;
    char operator!=(const SComplex& x) const {return 1 - ((*this) == x);}
};

SComplex operator+(const double& a, const SComplex& b);
SComplex operator-(const double& a, const SComplex& b);
SComplex operator*(const double& a, const SComplex& b);
SComplex operator/(const double& a, const SComplex& b);

SComplex operator==(const double& a, const SComplex& b);
SComplex operator!=(const double& a, const SComplex& b);

SComplex&makeZero(SComplex& x);
SComplex& makeOne(SComplex& x);
SComplex& makeInvalid(SComplex& x);

char isZero(const SComplex& x);
char isOne(const SComplex& x);

std::istream& operator>>(std::istream& input, SComplex& x);
std::ostream& operator<<(std::ostream& output, SComplex& x);

#endif
