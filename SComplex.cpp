#include "SComplex.h"

SComplex::SComplex(double x) {
    real = x;
    image = 0;
}

SComplex::SComplex(double a, double b) {
    real = a;
    image = b;
}

double SComplex::arg(void) const {
    double res;

    if ( isZero(this -> abs()) )
        res = 0;
    else if ( std::asin(this -> image / this -> abs()) > 0 )
        res = std::acos(this -> image / this -> abs());
    else
        res = -std::acos(this -> image / this -> abs());

    return res;
}

SComplex& SComplex::operator=(const double& x) {
    real = x;
    image = 0;

    return *this;
}

SComplex& SComplex::operator+=(const double& x) {
    real += x;

    return *this;
}

SComplex& SComplex::operator-=(const double& x) {
    real -= x;

    return *this;
}

SComplex& SComplex::operator*=(const double& x) {
    real *= x;
    image *= x;

    return *this;
}

SComplex& SComplex::operator/=(const double& x) {
    real /= x;
    image /= x;

    return *this;
}

SComplex SComplex::operator+(const double& x) const {
    SComplex res(*this);

    return res += x;
}

SComplex SComplex::operator-(const double& x) const {
    SComplex res(*this);

    return res -= x;
}

SComplex SComplex::operator*(const double& x) const {
    SComplex res(*this);

    return res *= x;
}

SComplex SComplex::operator/(const double& x) const {
    SComplex res(*this);

    return res /= x;
}

SComplex& SComplex::operator+=(const SComplex& x) {
    real += x.real;
    image += x.image;

    return *this;
}

SComplex& SComplex::operator-=(const SComplex& x) {
    real -= x.real;
    image -= x.image;

    return *this;
}

SComplex& SComplex::operator*=(const SComplex& x) {
    double tmp = real;

    real = real * x.real - image * x.image;
    image = x.image * tmp + x.real * image;

    return *this;
}

SComplex& SComplex::operator/=(const SComplex& x) {
    double tmp = real;
    double k = x.real * x.real + x.image * x.image;

    real = real * x.real + image * x.image;
    image = x.real * image - x.image * tmp ;

    return *this /= k;
}

SComplex SComplex::operator+(const SComplex& x) const {
    SComplex res(*this);

    return res += x;
}

SComplex SComplex::operator-(const SComplex& x) const {
    SComplex res(*this);

    return res -= x;
}

SComplex SComplex::operator*(const SComplex& x) const {
    SComplex res(*this);

    return res *= x;
}

SComplex SComplex::operator/(const SComplex& x) const {
    SComplex res(*this);

    return res /= x;
}

SComplex& SComplex::operator=(double& x) {
    real = x;
    image = 0;

    return *this;
}

SComplex& SComplex::operator~(void) {
    image *= -1;

    return *this;
}

char SComplex::operator==(const double& x) const {
    return cmp(real, x) * isZero(image);
}

char SComplex::operator==(const SComplex& x) const {
    return cmp(real, x.real) * cmp(image, x.image);
}

SComplex operator/(const double& a, const SComplex& b) {
    SComplex res(a);

    return res /= b;
}

SComplex operator+(const double& a, const SComplex& b) {
    return b + a;
}

SComplex operator-(const double& a, const SComplex& b) {
    return b - a;
}

SComplex operator*(const double& a, const SComplex& b) {
    return b * a;
}

SComplex operator==(const double& a, const SComplex& b) {
    return b == a;
}

SComplex operator!=(const double& a, const SComplex& b) {
    return b != a;
}

SComplex&makeZero(SComplex& x) {
    return x = 0;
}

SComplex& makeOne(SComplex& x) {
    return x = 1;
}

SComplex& makeInvalid(SComplex& x) {
    return x = NAN;
}

char isZero(const SComplex& x) {
    return isZero(x.abs());
}

char isOne(const SComplex& x) {
    return isOne(x.abs()) * isZero(x.arg());
}

std::istream& operator>>(std::istream& input, SComplex& x) {
    return input >> x.re() >> x.im();
}

std::ostream& operator<<(std::ostream& output, SComplex& x) {
    return output << x.re() << ' ' << x.im();
}
