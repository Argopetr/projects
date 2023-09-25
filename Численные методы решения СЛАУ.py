import random
import math

def vAbs(a, n):
    res = 0

    for i in range(n):
        if abs(a[i]) > res:
            res = abs(a[i])

    return res

# Reverse matrix to matrix A
# n - matrix dimension
# eps - precision
def reverse(M, n, eps):
    flag = 1
    copy = [[M[i][j] for j in range(n)] for i in range(n)]
    res = [[0 for j in range(n)] for i in range(n)]
    tmp1 = [0 for i in range(n)]
    tmp2 = [0 for i in range(n)]

    for i in range(n):
        res[i][i] = 1

    for i in range(n):
        flag = n
    
        while (abs(copy[i][i]) < eps and flag > 0):
            for k in range(n):
                tmp1[k] = copy[i][k]
                tmp2[k] = res[i][k]

            for j in range(i, n - 1):
                for k in range(n):
                    copy[j][k] = copy[j + 1][k]
                    res[j][k] = res[j + 1][k]

            for k in range(n):
                copy[n - 1][k] = tmp1[k]
                res[n - 1][k] = tmp2[k]

            flag -= 1

        while (abs(copy [i][i]) < 1):
            for j in range(n):
                copy [i][j] *= 2
                res[i][j] *= 2

        for k in range(i + 1, n):
            for j in range(n):
                res[k][j] -= res[i][j] * copy[k][i] / copy[i][i]

            for j in range(n - 1, i, -1):
                copy[k][j] -= copy[i][j] * copy[k][i] / copy[i][i]

    for i in range(n - 1, -1, -1):
        for k in range(i):
            for j in range(n):
                res[k][j] -= res[i][j] * copy[k][i] / copy[i][i]

        for j in range(n):
            res[i][j] /= copy[i][i]

    return res

# ||x2 - x1|| ^ 2
# n - vectors dimension
def iterationError(x2, x1, n):
    res = 0

    for i in range(n):
        for j in range(n):
            res += (x2[i] - x1[i]) * (x2[i] - x1[i])

    return res

# Solving of linear system by method of iteration
# A - matrix of system left part
# b - vector of system right part
# x0 - start solution
# w - iteration parameter
# eps - precision
def linearIterat(A, b, n, x0, w, eps):
    count = 0
    x2 = [0 for i in range(n)]

    if vAbs(b, n) > eps:
        D = reverse(A, n, eps)
        beta = [0 for i in range(n)]
        alfa = [[0 for j in range(n)] for i in range(n)]
        x1 = x0[:]
        x1[0] += 2 * eps

        count = 0

        for i in range(n):
            D[i][i] -= w
            x2[i] = x0[i]
            
            for j in range(n):
                beta[i] += D[i][j] * b[j]
                alfa[i][j] *= w

        while iterationError(x2, x1, n) >= eps * eps:
            x1 = x2[:]

            for i in range(n):
                x2[i] = beta[i]

                for j in range(n):
                    x2[i] += alfa[i][j] * x1[j]

            count += 1


    return [x2, count]

# ||A x0 - b|| ^ 2
# A - matrix from left part of system A x = b
# b - vector from right part
# n - matrix dimension
def solutionError(A, b, x0, n):
    res = 0
    temp = 0

    for i in range(n):
        temp = 0
        
        for j in range(n):
            temp += A[i][j] * x0[j]

        temp -= b[i]
        res += temp * temp

    return res

# Linear system solution by Successive over-relaxation method
# A - matrix of system left part
# b - vector of system right part
# N - system dimension
# eps - precision, w - relaxation parameter
def matrixRelax(A, b, n, eps, w):
    x1 = [0 for i in range(n)]

    count = 0

    if vAbs(b, n) > eps:
        x0 = [0 for i in range(n)]
        
        while solutionError(A, b, x1, n) >= eps * eps:
            x0 = x1[:]

            for i in range(n):
                x1[i] = (1 - w) * x0[i] + w * b[i] / A[i][i]

                for j in range(n):
                    if j < i:
                        x1[i] -= w * A[i][j] * x1[j] / A[i][i]
                    if j > i:
                        x1[i] -= w * A[i][j] * x0[j] / A[i][i]

            count += 1

    return [x1, count]

# (a, b)
# n - vectors dimension
def vMult(a, b, n):
    summ = 0

    for i in range(n):
        summ += a[i] * b[i]

    return summ

# Conjugate gradient method for linear system
# A - matrix, b - vector of system right part
# n - system dimension
# eps - precision
def conjGrad(A, b, n, eps): 
    x = [0 for i in range(n)]

    count = 0

    if vAbs(b, n) > eps:
        if vMult(b, b, n) > eps * eps: 
            r = [b[i] for i in range(n)]

            z = r[:]
            r0 = r[:]
            
            alfa = 1
            beta = 0
            
            while alfa * alfa * vMult(z, z, n) >= eps * eps:
                alfa = 0
                
                for i in range(n):
                    for j in range(n):
                        alfa += A[i][j] * z[i] * z[j]

                alfa = vMult(r0, r0, n) / alfa

                for i in range(n):
                    x[i] += alfa * z[i]

                    for j in range(n):
                        r[i] -= alfa * A[i][j] * z[j]

                beta = vMult(r, r, n) / vMult(r0, r0, n)
                r0 = r[:]

                for i in range(n):
                    z[i] = r[i] + beta * z[i]

                count += 1

    return [x, count]

# Biconjugate gradient method for linear system
# A - matrix, b - vector of system right part
# n - system dimension
# eps - precision
def biconjGrad(A, b, n, eps):
    x = [0 for i in range(n)]

    count = 0

    if vAbs(b, n) > eps:
        if vMult(b, b, n) > n * eps * eps:
            r0 = b[:]
            p0 = r0[:]
            p = r0[:]
            z = r0[:]
            s = r0[:]
            r = r0[:]

            alfa = 1
            beta = 0

            while alfa * alfa * vMult(z, z, n) >= eps * eps:
                alfa = vMult(p, z, n)
                beta = 0

                for i in range(n):
                    for j in range(n):
                        beta += A[i][j] * s[i] * z[j]

                alfa /= beta

                for i in range(n):
                    x[i] += alfa * z[i]

                    for j in range(n):
                        r[i] -= alfa * A[i][j] * z[j]
                        p[i] -= alfa * A[i][j] * s[j]

                beta = vMult(p, r, n) / vMult(p0, r0, n)

                for i in range(n):
                    z[i] = r[i] + beta * z[i]
                    s[i] = p[i] + beta * s[i]
                    r0[i] = r[i]
                    p0[i] = p[i]

                count += 1

    return [x, count]

# Lower triangular part of matrix
# N - matrix dimension
def L_LU(matrix, N):
    for k in range(N):
        if N - k > 0:
            for i in range(N - 1 - k):
                for j in range(N - 1 - k):
                    matrix[k + 1 + i][k + 1 + j] -= matrix[k + 1 + i][k] * matrix[k][k + 1 + j] / matrix[k][k]

            for i in range(N - 1 - k):
                matrix[k + 1 + i][k] /= matrix[k][k]
                matrix[k][k + 1 + i] = 0

        matrix[k][k] = 1

    return matrix

# Upper triangular part of matrix
# N - matrix dimension
def U_LU(matrix, N):
    for k in range(N):
        if N - k > 0:
            for i in range(N - 1 - k):
                for j in range(N - 1 - k):
                    matrix[k + 1 + i][k + 1 + j] -= matrix[k + 1 + i][k] * matrix[k][k + 1 + j] / matrix[k][k]

            for i in range(N - 1 - k):
                matrix[k + 1 + i][k] = 0

    return matrix

# Solution of linear system by LU-decomposition
# A - matrix of system left part
# b - vector of system right part
# N - system dimension
def LU_sol(A, b, N, eps):
    res = [0 for i in range(N)]

    if vAbs(b, N) > eps:
        L = [[0 for j in range(N)] for i in range(N)]
        U = [[0 for j in range(N)] for i in range(N)]

        for i in range(N):
            res[i] = b[i]
            
            for j in range(N):
                L[i][j] =  A[i][j]
                U[i][j] =  A[i][j]
        
        L_LU(L, N)
        U_LU(U, N)

        for i in range(N):
            for j in range(i):
                res[i] -= L[i][j] * res[j]

            res[i] /= L[i][i]

        for i in range(N):
            for j in range(i):
                res[N - 1 - i] -= U[N - 1 - i][N - 1 - j] *res[N - 1 - j]

            res[N - 1 - i] /= U[N - 1 - i][N - 1 - i]

    return res

eps = 0.000001
wIterat = 0.00001
wRelax = 0.5

wdth = 7
prc = 3

N = random.randint(1, 19) # Generation of space dimension
C = [[0 for j in range(N)] for i in range(N)] # Creation of transformation matrix
A = [[0 for j in range(N)] for i in range(N)] # Creation of matrix A from system Ax = b
D = [random.uniform(eps, 99) for i in range(N)] # Creation of eigenvalues vector
b = [random.uniform(-99, 99) for i in range(N)] # Generation of vector b from system Ax = b

tmp = 0

# Making orthogonal matrix C and positives eigenvalues from D
# A = CT D C, CT - transposed matrix C
for i in range(N):
    while vMult(C[i], C[i], N) < eps * eps:
        for j in range(N):
            C[i][j] = random.uniform(-99, 99)
        
        for j in range(N):
            tmp = vMult(C[i], C[j], N)

            for k in range(N):
                C[i][k] -= tmp * C[j][k]

    tmp = math.sqrt(vMult(C[i], C[i], N))

    for j in range(N):
        C[i][j] /= tmp

for i in range(N):
    for j in range(N):
        for k in range(N):
            A[i][j] += C[k][i] * D[k] * C[k][j]

s1 = linearIterat(A, b, N, b, wIterat, eps)
s2 = LU_sol(A, b, N, eps)
s3 = conjGrad(A, b, N, eps)
s4 = biconjGrad(A, b, N, eps)
s5 = matrixRelax(A, b, N, eps, wRelax)

print("Matrix A:")
for i in range(N):
    for j in range(N):
        print('{:{width}.{prec}f} '.format(A[i][j], width = wdth, prec = prc), end = '')
        
    print()

print()
print("vector b:")

for i in range(N):
    print('{:{width}.{prec}f}'.format(b[i], width = wdth, prec = prc))

print()
print("Solution of system Ax = b by method of iterations:")

for i in range(N):
    print('{:{width}.{prec}f}'.format(s1[0][i], width = wdth, prec = prc))  

print("Count of iterations:",  s1[1])
print()
print("Solution of system Ax = b by LU-decomposition method:")

for i in range(N):
    print('{:{width}.{prec}f}'.format(s2[i], width = wdth, prec = prc))

print()
print("Solution of system Ax = b by conjugate gradient method for linear system:")

for i in range(N):
    print('{:{width}.{prec}f}'.format(s3[0][i], width = wdth, prec = prc))

print("Count of iterations:",  s3[1])
print()
print("Solution of system Ax = b by biconjugate gradient method for linear system:")

for i in range(N):
    print('{:{width}.{prec}f}'.format(s4[0][i], width = wdth, prec = prc))

print("Count of iterations:",  s4[1])
print()
print("Solution of system Ax = b by successive over-relaxation method:")

for i in range(N):
    print('{:{width}.{prec}f}'.format(s5[0][i], width = wdth, prec = prc))
    
print("Count of iterations:",  s5[1])
