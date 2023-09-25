#ifndef LINAL_H
#define LINAL_H

#include <iostream>
#include "Neutral.h"

// VECTOR

template<typename Field>
struct SVector {
private:
    int n;
    Field *space;

public:
    SVector();
    SVector(int n);
    SVector(const SVector& v);
    ~SVector() {delete[] space;}

    SVector& zero(void);
    SVector& setDim(int i);

    int dim(void) const {return (n > 0)? n : -1;}
    char isValid(void) const {return (n > 0)? 1 : 0;}

    Field getItem(int i) const;
    char setItem(int i, const Field& val);

    Field& operator[](int i) {return space[i];}

    SVector& operator*=(const Field& k);
    SVector& operator/=(const Field& k);

    SVector operator*(const Field& k) const;
    SVector operator/(const Field& k) const;

    SVector& operator+=(const SVector& v);
    SVector& operator-=(const SVector& v);
    SVector& operator=(const SVector& v);

    SVector operator+(const SVector& v) const;
    SVector operator-(const SVector& v) const;

    Field operator*(const SVector& v) const;
};

template<typename Field>
SVector<Field>::SVector() {
        n = -1;
        space = NULL;
}

template<typename Field>
SVector<Field>::SVector(int n) {
    if ( n > 0 ) {
        this -> n = n;
        space = new Field[n];

        for(int i; i < n; i++)
            makeZero(space[i]);
    } else {
        n = -1;
        space = NULL;
    }
}

template<typename Field>
SVector<Field>::SVector(const SVector<Field>& v) {
    if ( &v != this ) {
        if ( v.n > 0 ) {
            int i = v.n;

            n = v.n;
            space = new Field[n];

            for(i--; i >= 0; i--)
                space[i] = v.space[i];
        } else {
            n = -1;
            space = NULL;
        }
    }
}

template<typename Field>
SVector<Field>& SVector<Field>::zero(void) {
    for(int i = 0; i < n; i++)
        makeZero(space[i]);

    return *this;
}

template<typename Field>
SVector<Field>& SVector<Field>::setDim(int i) {
    Field* newSpace = NULL;

    if ( i > 0 ) {
        newSpace = new Field[i];

        for(int k = 0; k < i; k++)
            if ( k < n )
                newSpace[k] = space[k];
            else
                makeZero(newSpace[k]);
    } else
        i = -1;

    delete[] space;

    n = i;
    space = newSpace;

    return *this;
}

template<typename Field>
Field SVector<Field>::getItem(int i) const {
    Field res;

    if ( i >= 0 && i < n )
        res = space[i];
    else
        makeInvalid(res);

    return res;
}

template<typename Field>
char SVector<Field>::setItem(int i, const Field& val) {
    char err = 0;

    if ( i >= 0 && i < n )
        space[i] = val;
    else
        err = -1;

    return err;
}

template<typename Field>
SVector<Field>& SVector<Field>::operator*=(const Field& k) {
    int i = n;

    for(i = 0; i < n; i++)
        space[i] *= k;

    return *this;
}

template<typename Field>
SVector<Field>& SVector<Field>::operator/=(const Field& k) {
    int i = n;

    for(i = 0; i < n; i++)
        space[i] /= k;

    return *this;
}

template<typename Field>
SVector<Field> SVector<Field>::operator*(const Field& k) const {
    SVector<Field> res(*this);

    return res *= k;
}

template<typename Field>
SVector<Field> SVector<Field>::operator/(const Field& k) const {
    SVector<Field> res(*this);

    return res /= k;
}

template<typename Field>
SVector<Field>& SVector<Field>::operator=(const SVector& v) {
    if ( &v != this ) {
        delete[] space;

        if ( v.n > 0 ) {
            int i = v.n;

            n = v.n;
            space = new Field[n];

            for(i--; i >= 0; i--)
                space[i] = v.space[i];
        } else {
            n = -1;
            space = NULL;
        }
    }

    return *this;
}

template<typename Field>
SVector<Field>& SVector<Field>::operator+=(const SVector<Field>& v) {
    if ( v.n == n )
        for(int i = 0; i < n; i++)
            space[i] += v.space[i];
    else {
        delete[] space;

        n = -1;
        space = NULL;
    }

    return *this;
}

template<typename Field>
SVector<Field>& SVector<Field>::operator-=(const SVector<Field>& v) {
    if ( v.n == n )
        for(int i = 0; i < n; i++)
            space[i] -= v.space[i];
    else {
        delete[] space;

        n = -1;
        space = NULL;
    }

    return *this;
}

template<typename Field>
SVector<Field> SVector<Field>::operator+(const SVector<Field>& v) const {
    SVector<Field> res(*this);

    return res += v;
}

template<typename Field>
SVector<Field> SVector<Field>::operator-(const SVector<Field>& v) const {
    SVector<Field> res(*this);

    return res -= v;
}

template<typename Field>
Field SVector<Field>::operator*(const SVector<Field>& v) const {
    Field res;

    if ( v.n == n ) {
        makeZero(res);

        for(int i = 0; i < n; i++)
            res += space[i] * v.space[i];
    } else
        makeInvalid(res);

    return res;
}

template<typename Field>
SVector<Field> operator*(const Field& k, const SVector<Field>& v) {return v * k;}

// MATRIX

template<typename Field>
struct SMatrix {
private:
    int n, m;
    SVector<Field> *space;

public:
    SMatrix();
    SMatrix(int i);
    SMatrix(int i, int j);
    SMatrix(const SMatrix& M);
    ~SMatrix() {delete[] space;}

    SMatrix& E(void);
    SMatrix& zero(void);
    SMatrix& setDim(int i, int j);

    Field det(void) const;
    int rank(void) const;
    int stringDim(void) const {return (n > 0 && m > 0)? n : -1;}
    int columnDim(void) const {return (n > 0 && m > 0)? m : -1;}
    char isValid(void) const {return (n > 0 && m > 0)? 1 : 0;}

    char setItem(int i, int j, Field val);
    Field getItem(int i, int j) const;

    char setString(int i, const SVector<Field>& v);
    char setColumn(int j, const SVector<Field>& v);

    SVector<Field> getString(int i) const;
    SVector<Field> getColumn(int j) const;

    SMatrix trans(void) const;
    SMatrix inv(void) const;

    SMatrix nullSpace(void) const;
    SVector<Field> solve(const SVector<Field>& v) const;

    SVector<Field>& operator[](int i) {return space[i];}

    SMatrix& operator*=(const Field& k);
    SMatrix& operator/=(const Field& k);

    SMatrix operator*(const Field& k) const;
    SMatrix operator/(const Field& k) const;

    SMatrix& operator=(const SMatrix& M);
    SMatrix& operator+=(const SMatrix& M);
    SMatrix& operator-=(const SMatrix& M);
    SMatrix& operator*=(const SMatrix& M);

    SMatrix operator+(const SMatrix& M) const;
    SMatrix operator-(const SMatrix& M) const;
    SMatrix operator*(const SMatrix& M) const;

    SVector<Field> operator*(const SVector<Field>& v) const;
};

template<typename Field>
SMatrix<Field>::SMatrix() {
    n = -1;
    m = -1;
    space = NULL;
}

template<typename Field>
SMatrix<Field>::SMatrix(int i) {
    if ( i > 0 ) {

        n = i;
        m = i;

        space = new SVector<Field>[i];

        for(int k = 0; k < i; k++) {
            space[k].setDim(i);
            makeOne(space[k][k]);
        }
    } else {
        n = -1;
        m = -1;
        space = NULL;
    }
}

template<typename Field>
SMatrix<Field>::SMatrix(int i, int j) {
    if ( i > 0 && j > 0 ) {
        SVector<Field> tmp(j);

        n = i;
        m = j;

        space = new SVector<Field>[i];

        for(int k = 0; k < i; k++)
            space[k] = tmp;
    } else {
        n = -1;
        m = -1;
        space = NULL;
    }
}

template<typename Field>
SMatrix<Field>::SMatrix(const SMatrix<Field>& M) {
    if ( M.n > 0 && M.m > 0 ) {
        n = M.n;
        m = M.m;

        space = new SVector<Field>[n];

        for(int k = 0; k < n; k++)
            space[k] = M.space[k];
    } else {
        n = -1;
        m = -1;
        space = NULL;
    }
}

template<typename Field>
SMatrix<Field>& SMatrix<Field>::E(void) {
    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            if ( i == j )
                makeOne(space[i][j]);
            else
                makeZero(space[i][j]);

    return *this;
}

template<typename Field>
SMatrix<Field>& SMatrix<Field>::zero(void) {
    for(int i = 0; i < n; i++)
        space[i].zero();

    return *this;
}

template<typename Field>
SMatrix<Field>& SMatrix<Field>::setDim(int i, int j) {
    SVector<Field> *newSpace = NULL;

    if ( i > 0 && j > 0 ) {
        newSpace = new SVector<Field>[i];

        for(int k = 0; k < i; k++)
            if ( k < n )
                (newSpace[k] = space[k]).setDim(j);
            else
                newSpace[k].setDim(j);
    } else {
        i = -1;
        j = -1;
    }

    delete[] space;

    space = newSpace;
    n = i;
    m = j;

    return *this;
}

template<typename Field>
Field SMatrix<Field>::det(void) const {
    Field res;

    if ( n > 0 && m > 0 && n == m ) {
        int i, j, k;
        char flag = 1;
        SMatrix<Field> copy(*this);

        for(i = 0; i < n && flag; i++)
            if ( isZero(copy.space[i][i]) )
                flag = 0;
            else
                for(k = i + 1; k < n; k++)
                    for(j = m - 1; j >= i; j--)
                        copy.space[k][j] -= copy.space[i][j] * copy.space[k][i] / copy.space[i][i];

        if ( flag )
            for(makeOne(res), k = 0; k < n; k++)
                res *= copy.space[k][k];
        else
            makeZero(res);

    } else
        makeInvalid(res);

    return res;
}

template<typename Field>
int SMatrix<Field>::rank(void) const {
    int res;

    if ( n > 0 && m > 0 ) {
        int i, j, k, flag, shift = 0;
        SVector<Field> tmp;
        SMatrix<Field> copy(*this);

        res = (n < m)? n : m;

        for(i = 0; i < (n < m)? n : m && res; i++) {
            flag = n - i;

            for(; flag && isZero(copy.space[i][i + shift]); flag--) {
                tmp = copy.space[i];

                for(k = i; k < n - 1; k++)
                    copy.space[k] = copy.space[k + 1];

                copy.space[n - 1] = tmp;
            }

            if ( flag )
                for(k = i + 1; k < n; k++)
                    for(j = m - 1; j >= i + shift; j--)
                        copy.space[k][j] -= copy.space[i][j] * copy.space[k][i + shift] / copy.space[i][i + shift];
            else {
                shift++;
                res--;
            }
        }
    } else
        res = -1;

    return res;
}

template<typename Field>
char SMatrix<Field>::setItem(int i, int j, Field val) {
    char err = 0;

    if ( n > 0 && m > 0 && i >= 0 && j >= 0 && i < n && j < m )
        space[i][j] = val;
    else
        err = -1;

    return err;
}

template<typename Field>
Field SMatrix<Field>::getItem(int i, int j) const {
    Field res;

    if ( n > 0 && m > 0 && i >= 0 && j >= 0 && i < n && j < m )
        res = space[i].getItem(j);
    else
        makeInvalid(res);

    return res;
}

template<typename Field>
char SMatrix<Field>::setString(int i, const SVector<Field>& v) {
    char err = 0;

    if ( n > 0 && m > 0 && v.dim() == m && i >= 0 && i < n )
        space[i] = v;
    else
        err = -1;

    return err;
}

template<typename Field>
char SMatrix<Field>::setColumn(int j, const SVector<Field>& v) {
    char err = 0;

    if ( n > 0 && m > 0 && j >= 0 && j < 0 && v.dim() == n )
        for(int i = 0; i < n; i++)
            space[i][j] = v.getItem(i);
    else
        err = -1;

    return err;
}

template<typename Field>
SVector<Field> SMatrix<Field>::getString(int i) const {
    SVector<Field> res;

    if ( n > 0 && m > 0 && i >= 0 && i < n )
        res = space[i];

    return res;
}

template<typename Field>
SVector<Field> SMatrix<Field>::getColumn(int j) const {
    SVector<Field> res(m * ((n > 0 && m > 0 && j >= 0 && j < m)? 1 : 0));

    if ( n > 0 && m > 0 && j >= 0 && j < m )
        for(int i = 0; i < n; i++)
            res[i] = space[i].getItem(j);

    return res;
}

template<typename Field>
SMatrix<Field> SMatrix<Field>::trans(void) const {
    SMatrix<Field> res(m, n);

    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            res.space[j][i] = space[i].getItem(j);

    return res;
}

template<typename Field>
SMatrix<Field> SMatrix<Field>::inv(void) const {
    SMatrix<Field> res;

    if ( n > 0 && m > 0 && n == m ) {
        int i, j, k, flag = 1;
        SMatrix<Field> tmp(n, n);
        SMatrix<Field> copy(*this);

        tmp.E();

        for(i = 0; i < n && flag; i++)
            if ( isZero(copy.space[i][i]) )
                flag = 0;
            else
                for(k = n - 1; k > i; k--) {
                    for(j = m - 1; j >= 0; j--)
                        tmp.space[k][j] -= tmp.space[i][j] * copy.space[k][i] / copy.space[i][i];

                    for(j = m - 1; j >= i; j--)
                        copy.space[k][j] -= copy.space[i][j] * copy.space[k][i] / copy.space[i][i];
                }

        if ( flag ) {
            for(i = n - 1; i >= 0; i--) {
                for(k = i - 1; k >= 0; k--)
                    for(j = m - 1; j >= 0; j--)
                        tmp.space[k][j] -= tmp.space[i][j] * copy.space[k][i] / copy.space[i][i];

                for(j = m - 1; j >= 0; j--)
                    tmp.space[i][j] /= copy.space[i][i];
            }

            res -> n = res -> m = n;
            res -> space = tmp.space;
            tmp.space = NULL;
        }
    }

    return res;
}

template<typename Field>
SMatrix<Field> SMatrix<Field>::nullSpace(void) const {
    SMatrix<Field> res;

    if ( n > 0 && m > 0 ) {
        int j, k, flag, dim, rk;
        int *generalVar = new int[n];
        int *freeVar = new int[m];
        SVector<Field> tmp;
        SMatrix<Field> Mcopy(*this);

        for(dim = 0, rk = 0; rk + dim < m && rk < n;) {
            flag = n - rk;

            for(; flag && isZero(Mcopy.space[rk][rk + dim]); flag--) {
                tmp = Mcopy.space[rk];

                for(k = rk; k < n - 1; k++)
                    Mcopy.space[k] = Mcopy.space[k + 1];

                Mcopy.space[n - 1] = tmp;
            }

            if ( flag ) {
                for(k = 0; k < n; k++)
                    if ( k != rk )
                        for(j = m - 1; j >= rk + dim; j--)
                            Mcopy.space[k][j] -= Mcopy.space[rk][j] * Mcopy.space[k][rk + dim] / Mcopy.space[rk][rk + dim];

                for(j = m - 1; j >= rk + dim; j--)
                    Mcopy.space[rk][j] /= Mcopy.space[rk][rk + dim];

                generalVar[rk] = rk + dim;
                rk++;
            } else {
                freeVar[dim] = rk + dim;
                dim++;
            }
        }

        for(; rk + dim < m; dim++)
            freeVar[dim] = rk + dim;

        res.setDim(m, dim);

        for(j = 0; j < dim; j++) {
            makeOne(res.space[freeVar[j]][j]);

            for(k = 0; k < rk; k++)
                    res.space[generalVar[k]][j] -= Mcopy.space[k][freeVar[j]];
        }

	delete[] generalVar;
	delete[] freeVar;
    }

    return res;
}

template<typename Field>
SVector<Field> SMatrix<Field>::solve(const SVector<Field>& v) const {
    SVector<Field> res;

    if ( n > 0 && m > 0 && v.dim() == n ) {
        int i, j, k, flag, dim, rk;
        SMatrix<Field> Mcopy(*this);
        SVector<Field> Vcopy(v), tmp;
        Field bff;

        for(rk = 0, dim = 0; rk + dim < m && rk < n;) {
            flag = n - rk;

            for(; flag && isZero(Mcopy.space[rk][rk + dim]); flag--) {
                tmp = Mcopy.space[rk];
                bff = Vcopy[rk];

                for(k = rk; k < n - 1; k++){
                    Vcopy[k] = Vcopy[k + 1];
                    Mcopy.space[k] = Mcopy.space[k + 1];
                }

                Vcopy[n - 1] = bff;
                Mcopy.space[n - 1] = tmp;
            }

            if ( flag ) {
                for(k = 0; k < n; k++)
                    if ( k != rk ) {
                        Vcopy[k] -= Vcopy[rk] * Mcopy.space[k][rk + dim] / Mcopy.space[rk][rk + dim];

                        for(j = m - 1; j >= rk + dim; j--)
                            Mcopy.space[k][j] -= Mcopy.space[rk][j] * Mcopy.space[k][rk + dim] / Mcopy.space[rk][rk + dim];
                    }

                rk++;
            } else
                dim++;
        }

        flag = 1;

        for(i = rk; i < n && flag; i++)
            if ( !isZero(Vcopy[i]) )
                flag = 0;

        if ( flag ) {
            res.setDim(m);

            for(i = 0; i < rk; i++) {
                for(j = i; isZero(Mcopy.space[i][j]); j++);

                res[j] = Vcopy[i] / Mcopy.space[i][j];
            }
        }
    }

    return res;
}

template<typename Field>
SMatrix<Field>& SMatrix<Field>::operator*=(const Field& k) {
    for(int i = 0; i < n; i++)
        space[i] *= k;

    return *this;
}

template<typename Field>
SMatrix<Field>& SMatrix<Field>::operator/=(const Field& k) {
    for(int i = 0; i < n; i++)
        space[i] /= k;

    return *this;
}

template<typename Field>
SMatrix<Field> SMatrix<Field>::operator*(const Field& k) const {
    SMatrix<Field> res(*this);

    return res *= k;
}

template<typename Field>
SMatrix<Field> SMatrix<Field>::operator/(const Field& k) const {
    SMatrix<Field> res(*this);

    return res /= k;
}

template<typename Field>
SMatrix<Field>& SMatrix<Field>::operator=(const SMatrix<Field>& M) {
    if ( &M != this ) {
        delete[] space;

        if ( n > 0 && m > 0 ) {
            n = n;
            m = m;

            space = new SVector<Field>[n];

            for(int k = 0; k < n; k++)
                space[k] = M.space[k];
        } else {
            n = -1;
            m = -1;
            space = NULL;
        }
    }

    return *this;
}

template<typename Field>
SMatrix<Field>& SMatrix<Field>::operator+=(const SMatrix<Field>& M) {
    if ( n > 0 && m > 0 && n == n && m == m )
        for(int i = 0; i < n; i++)
            space[i] += M.space[i];
    else {
        delete[] space;

        n = -1;
        m = -1;
        space = NULL;
    }

    return *this;
}

template<typename Field>
SMatrix<Field>& SMatrix<Field>::operator-=(const SMatrix<Field>& M) {
    if ( n > 0 && m > 0 && n == n && m == m )
        for(int i = 0; i < n; i++)
            space[i] -= M.space[i];
    else {
        delete[] space;

        n = -1;
        m = -1;
        space = NULL;
    }

    return *this;
}

template<typename Field>
SMatrix<Field>& SMatrix<Field>::operator*=(const SMatrix<Field>& M) {
    if ( n > 0 && m > 0 && m == n && m > 0 ) {
        SMatrix<Field> tmp(n, m);

        for(int i = 0; i < n; i++)
            for(int j = 0; j < m; j++)
                for(int k = 0; k < m; k++)
                    tmp.space[i][j] += space[i].getItem(k) * M.space[k].getItem(j);

        delete[] space;

        space = tmp.space;
        tmp.space = NULL;
        m = m;
    } else {
        delete[] space;

        n = -1;
        m = -1;
    }

    return *this;
}

template<typename Field>
SMatrix<Field> SMatrix<Field>::operator+(const SMatrix<Field>& M) const {
    SMatrix<Field> res(*this);

    return res += M;
}

template<typename Field>
SMatrix<Field> SMatrix<Field>::operator-(const SMatrix<Field>& M) const {
    SMatrix<Field> res(*this);

    return res -= M;
}

template<typename Field>
SMatrix<Field> SMatrix<Field>::operator*(const SMatrix<Field>& M) const {
    SMatrix<Field> res(*this);

    return res *= M;
}

template<typename Field>
SVector<Field> SMatrix<Field>::operator*(const SVector<Field>& v) const {
    SVector<Field> res(((n > 0 && m > 0 && v.dim() == m)? 1 : 0) * n);

    for(int i = 0; i < n; i++)
        res.setItem(i, space[i] * v);

    return res;
}

template<typename Field>
SMatrix<Field> operator*(const Field& k, const SMatrix<Field>& M) {return M * k;}

template<typename Field>
SVector<Field>& operator*=(SVector<Field>& v, const SMatrix<Field>& M) {
    if ( M.columnDim() > 0 && M.stringDim() == v.dim() && v.dim() > 0 ) {
        SVector<Field> tmp(M.columnDim());

        for(int i = 0; i < v.dim(); i++)
            for(int j = 0; j < M.columnDim(); j++)
                tmp[j] += v[i] * M.getItem(i, j);

        v = tmp;
    } else
        v.setDim(-1);

    return v;
}

template<typename Field>
SVector<Field> operator*(SVector<Field>& v, const SMatrix<Field>& M) {
    SVector<Field> res(v);

    return res *= M;
}

template<typename Field>
std::istream& operator>>(std::istream& input, SVector<Field>& v) {
    int n;
    Field tmp;

    input >> n;
    v.setDim(n);

    for(int i = 0; i < n; i++) {
        input >> tmp;
        v[i] = tmp;
    }

    return input;
}

template<typename Field>
std::istream& operator>>(std::istream& input, SMatrix<Field>& M) {
    int n, m;
    Field tmp;

    input >> n >> m;
    M.setDim(n, m);

    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++) {
            input >> tmp;
            M.setItem(i, j, tmp);
        }

    return input;
}

template<typename Field>
std::ostream& operator<<(std::ostream& output, const SVector<Field>& v) {
    if ( v.isValid() ) {
        Field tmp;
        int i;

        for(i = 0; i < v.dim(); i++) {
            tmp = v.getItem(i);
            output << tmp << ' ';
        }
    } else
        output << "NAV";

    return output;
}

template<typename Field>
std::ostream& operator<<(std::ostream& output, const SMatrix<Field>& M) {
    if ( M.isValid() ) {
        int i, j;
        Field tmp;

        for(i = 0; i < M.stringDim(); i++) {
            for(j = 0; j < M.columnDim(); j++) {
                tmp = M.getItem(i, j);
                output << tmp << ' ';
            }

            output << std::endl;
        }
    } else
        output << "NAM";

    return output;
}

#endif
