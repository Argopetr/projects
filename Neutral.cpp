#include "Neutral.h"

#define NEUTRAL_EPS 0.000001

#define MAX(a, b) ((a > b)? (a) : (b))
#define IS_ZERO(x) (std::fabs(x) < NEUTRAL_EPS)
#define CMP(a, b) ((std::fabs((a) - (b)) < NEUTRAL_EPS * MAX(std::fabs(a), std::fabs(b)))? 1 : 0)

char cmp(long double x, long double y) {return CMP(x, y);}

// CHAR

char& makeZero(char& x) {return x = 0;}
char& makeOne(char& x) {return x = 1;}

char isZero(const char& x) {return (x == 0)? 1 : 0;}
char isOne(const char& x) {return (x == 1)? 1 : 0;}

// UNSIGNED CHAR

unsigned char& makeZero(unsigned char& x) {return x = 0;}
unsigned char& makeOne(unsigned char& x) {return x = 1;}

char isZero(const unsigned char& x) {return (x == 0)? 1 : 0;}
char isOne(const unsigned char& x) {return (x == 1)? 1 : 0;}

// SIGNED CHAR

signed char& makeZero(signed char& x) {return x = 0;}
signed char& makeOne(signed char& x) {return x = 1;}

char isZero(const signed char& x) {return (x == 0)? 1 : 0;}
char isOne(const signed char& x) {return (x == 1)? 1 : 0;}

// SHORT

short& makeZero(short& x) {return x = 0;}
short& makeOne(short& x) {return x = 1;}

char isZero(const short& x) {return (x == 0)? 1 : 0;}
char isOne(const short& x) {return (x == 1)? 1 : 0;}

// UNSIGNED SHORT

unsigned short& makeZero(unsigned short& x) {return x = 0;}
unsigned short& makeOne(unsigned short& x) {return x = 1;}

char isZero(const unsigned short& x) {return (x == 0)? 1 : 0;}
char isOne(const unsigned short& x) {return (x == 1)? 1 : 0;}

// INT

int& makeZero(int& x) {return x = 0;}
int& makeOne(int& x) {return x = 1;}

char isZero(const int& x) {return (x == 0)? 1 : 0;}
char isOne(const int& x) {return (x == 1)? 1 : 0;}

// UNSIGNED INT

unsigned int& makeZero(unsigned int& x) {return x = 0;}
unsigned int& makeOne(unsigned int& x) {return x = 1;}

char isZero(const unsigned int& x) {return (x == 0)? 1 : 0;}
char isOne(const unsigned int& x) {return (x == 1)? 1 : 0;}

// LONG INT

long int& makeZero(long int& x) {return x = 0;}
long int& makeOne(long int& x) {return x = 1;}

char isZero(const long int& x) {return (x == 0)? 1 : 0;}
char isOne(const long int& x) {return (x == 1)? 1 : 0;}

// UNSIGNED LONG INT

unsigned long int& makeZero(unsigned long int& x) {return x = 0;}
unsigned long int& makeOne(unsigned long int& x) {return x = 1;}

char isZero(const unsigned long int& x) {return (x == 0)? 1 : 0;}
char isOne(const unsigned long int& x) {return (x == 1)? 1 : 0;}

// LONG LONG INT

long long int& makeZero(long long int& x) {return x = 0;}
long long int& makeOne(long long int& x) {return x = 1;}

char isZero(const long long int& x) {return (x == 0)? 1 : 0;}
char isOne(const long long int& x) {return (x == 1)? 1 : 0;}

// UNSIGNED LONG LONG INT

unsigned long long int& makeZero(unsigned long long int& x) {return x = 0;}
unsigned long long int& makeOne(unsigned long long int& x) {return x = 1;}

char isZero(const unsigned long long int& x) {return (x == 0)? 1 : 0;}
char isOne(const unsigned long long int& x) {return (x == 1)? 1 : 0;}

// FLOAT

float& makeZero(float& x) {return x = 0;}
float& makeOne(float& x) {return x = 1;}
float& makeInvalid(float& x) {return x = NAN;}

char isZero(const float& x) {return IS_ZERO(x)? 1 : 0;}
char isOne(const float& x) {return CMP(x, 1)? 1 : 0;}

// DOUBLE

double& makeZero(double& x) {return x = 0;}
double& makeOne(double& x) {return x = 1;}
double& makeInvalid(double& x) {return x = NAN;}

char isZero(const double& x) {return IS_ZERO(x)? 1 : 0;}
char isOne(const double& x) {return CMP(x, 1)? 1 : 0;}

// LONG DOUBLE

long double& makeZero(long double& x) {return x = 0;}
long double& makeOne(long double& x) {return x = 1;}
long double& makeInvalid(long double& x) {return x = NAN;}

char isZero(const long double& x) {return IS_ZERO(x)? 1 : 0;}
char isOne(const long double& x) {return CMP(x, 1)? 1 : 0;}
