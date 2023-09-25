#ifndef NEUTRAL_H
#define NEUTRAL_H

#include <cmath>

char cmp(long double x, long double y);

// CHAR

char& makeZero(char& x);
char& makeOne(char& x);

char isZero(const char& x);
char isOne(const char& x);

// UNSIGNED CHAR

unsigned char& makeZero(unsigned char& x);
unsigned char& makeOne(unsigned char& x);

char isZero(const unsigned char& x);
char isOne(const unsigned char& x);

// SIGNED CHAR

signed char& makeZero(signed char& x);
signed char& makeOne(signed char& x);

char isZero(const signed char& x);
char isOne(const signed char& x);

// SHORT

short& makeZero(short& x);
short& makeOne(short& x);

char isZero(const short& x);
char isOne(const short& x);

// UNSIGNED SHORT

unsigned short& makeZero(unsigned short& x);
unsigned short& makeOne(unsigned short& x);

char isZero(const unsigned short& x);
char isOne(const unsigned short& x);

// INT

int& makeZero(int& x);
int& makeOne(int& x);

char isZero(const int& x);
char isOne(const int& x);

// UNSIGNED INT

unsigned int& makeZero(unsigned int& x);
unsigned int& makeOne(unsigned int& x);

char isZero(const unsigned int& x);
char isOne(const unsigned int& x);

// LONG INT

long int& makeZero(long int& x);
long int& makeOne(long int& x);

char isZero(const long int& x);
char isOne(const long int& x);

// UNSIGNED LONG INT

unsigned long int& makeZero(unsigned long int& x);
unsigned long int& makeOne(unsigned long int& x);

char isZero(const unsigned long int& x);
char isOne(const unsigned long int& x);

// LONG LONG INT

long long int& makeZero(long long int& x);
long long int& makeOne(long long int& x);

char isZero(const long long int& x);
char isOne(const long long int& x);

// UNSIGNED LONG LONG INT

unsigned long long int& makeZero(unsigned long long int& x);
unsigned long long int& makeOne(unsigned long long int& x);

char isZero(const unsigned long long int& x);
char isOne(const unsigned long long int& x);

// FLOAT

float& makeZero(float& x);
float& makeOne(float& x);
float& makeInvalid(float& x);

char isZero(const float& x);
char isOne(const float& x);

// DOUBLE

double& makeZero(double& x);
double& makeOne(double& x);
double& makeInvalid(double& x);

char isZero(const double& x);
char isOne(const double& x);

// LONG DOUBLE

long double& makeZero(long double& x);
long double& makeOne(long double& x);
long double& makeInvalid(long double& x);

char isZero(const long double& x);
char isOne(const long double& x);

#endif
