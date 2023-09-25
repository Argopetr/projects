import math

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
    D = reverse(A, n, eps)
    beta = [0 for i in range(n)]
    alfa = [[0 for j in range(n)] for i in range(n)]
    x2 = x0[:]
    x1 = x2[:]
    x1[0] += 2 * eps

    for i in range(n):
        D[i][i] -= w
        
        for j in range(n):
            beta[i] += D[i][j] * b[j]
            alfa[i][j] *= w

    while iterationError(x2, x1, n) >= eps * eps:
        x1 = x2[:]

        for i in range(n):
            x2[i] = beta[i]

            for j in range(n):
                x2[i] += alfa[i][j] * x1[j]

    return x2

# x0 + dx
# x0 - vector of start argument
# n - dimention
# dx - difference between new and start arguments
# i - index of difference direction
def step(x0, n, dx, i):
    res = x0[:]
    res[i] += dx
    
    return res

# A - B
# n - dimention of matrixes A, B
def matrixDiff(A, B, n):
    return [[A[i][j] - B[i][j] for j in range(n)] for i in range(n)]

# Biggest abs of matrix elements
# n - dimension of A matrix
def matrixAbs(A, n):
    res = 0

    for i in range(n):
        for j in range(n):
            if abs(A[i][j]) > res:
                res = abs(A[i][j])

    return res

# Biggest abs of vector elements
# n - dimension of a vector
def vectorAbs(a, n):
    res = 0

    for i in range(n):
        if abs(a[i]) > res:
            res = abs(a[i])

    return res

# (a, a)
# n - vector a dimension
def vLenSqr(a, n):
    summ = 0

    for i in range(n):
        summ += a[i] * a[i]

    return summ

# Jacobi matrix of function f(x), x = x0
# n - function dimension
# h0 - start step
# e - precision
def jacobi(f, n, x0, h0, eps):
    x = [x0[:] for i in range(n)]
    h = h0

    W0 = [[0 for j in range(n)] for i in range(n)]
    W1 = [[(f(step(x0, n, h, j))[i] - f(x0)[i]) / h for j in range(n)] for i in range(n)]

    while matrixAbs(matrixDiff(W1, W0, n), n) >= eps:
        h /= 2
        
        for i in range(n):
            for j in range(n):
                W0[i][j] = W1[i][j]
                W1[i][j] = (f(step(x0, n, h, j))[i] - f(x0)[i]) / h

    return W1

# Solution of non-linear system by secant method
# f - vector of function from system f(x) = 0
# x0 - start solution
# n - function dimension
# e - precision
def secantSol(f, n, x0, w, eps):
    A0 = [[0 for j in range(n)] for i in range(n)]
    A = jacobi(f, n, x0, 1, eps)
    s = linearIterat(A, f(x0), n, f(x0), w, eps)
    sTs = vLenSqr(s, n)
    x1 = x0[:]
    x2 = x0[:]

    for i in range(n):
        s[i] *= -1
        x2[i] += s[i]

    while sTs >= eps * eps:
        for i in range(n):
            for j in range(n):
                A0[i][j] = A[i][j]

        for i in range(n):
            for j in range(n):
                A[i][j] = A0[i][j] + (f(x2)[i] - f(x1)[i]) * s[j] / sTs

                for k in range(n):
                    A[i][j] -= A0[i][k] * s[k] * s[j] / sTs

        s = linearIterat(A, f(x2), n, f(x2), w, eps)
        sTs = vLenSqr(s, n)

        for i in range(n):
            x1[i] = x2[i]
            s[i] *= -1
            x2[i] += s[i]

    return x2

# Solution of non-linear system by Newton's method
# f - vector of function from system f(x) = 0
# x0 - start solution
# n - function dimension
# e - precision
def NewtonSol(f, n, x0, w, eps):
    x = x0[:]
    s = [0 for i in range(n)]
    s[0] = 2 * eps

    while vectorAbs(s, n) >= eps * eps:
        s = linearIterat(jacobi(f, n, x, 1, eps), f(x), n, f(x), w, eps)

        for i in range(n):
            s[i] *= -1
            x[i] += s[i]
    
    return x

# Runge-Kutta methods of Cauchy problem
# n - dimension, f - vector of functions from system u' = f(x, u)
# x0, u0 - initial condition u(x0) = u0
# x - argument of shearched value u(x)
# a, b, c - Butcher tableau
# s - number of stages, h - step of iteration
def rungeKutta(n, f, x0, u0, x, a, b, c, s, h):
    k = [[0 for j in range(n)] for i in range(s)]
    v_arg = [0 for i in range(s)]
    y1 = u0[:]
    y2 = u0[:]

    xi = x0
    delta = h
    m = 0
    
    while m * h < abs(x - x0):
        m += 1

    if x - x0 < 0:
        delta *= -1

    for i in range(m):
        if i >= m - 2:
            delta = x - xi
        
        for p in range(s):
            v_arg = y1[:]
                
            for t in range(n):
                for q in range(p):
                    v_arg[t] += delta * a[p - 1][q] * k[q][t]
                    
            k[p] = f(xi + c[p] * delta, v_arg)

            for t in range(n):
                y2[t] +=  delta * b[p] * k[p][t]
        
        xi += delta
        y1 = y2[:]

    return y2

# Runge-Kutta methods (s = 4)
def rungeKutta4(n, f, x0, u0, x, h):
    c = [0, 0.5, 0.5, 1]
    a = [[0.5], [0, 0.5], [0, 0, 1]]
    b = [1 / 6, 1 / 3, 1 / 3, 1 / 6]

    return rungeKutta(n, f, x0, u0, x, a, b, c, 4, h)

# Solving of Cauchy's problem with parameters
# F - vector of functions from system u' = F(x, u)
# n - dimension of F
# a - start argument
# alfa - vectror of functions (alfa0, alfa1, ..., alfa(k - 1))
# If u(a) = u0, then u0 = (alfa0(p), alfa1(p), ..., alfa(k - 1)(p), p0, p1, ..., p(n-1))
# Dimension of alfa is k
# p - vector of last n - k components of vector u0
# Dimension of p is n - k
# x - argument of searched value u(x)
# h - step of iteration
def CauchyProblemParameters(F, n, a, alfa, k, p, x, h):
    u0 = [0 for i in range(n)]

    for i in range(n):
        if i < k:
            u0[i] = alfa(p)[i]
        else:
            u0[i] = p[i - k]

    return rungeKutta4(n, F, a, u0, x, h)

# Solving of border problem
# F - vector of functions from system u' = F(x, u)
# n - dimension of F
# a - start argument
# b - argumet of stitching
# alfa - vectror of functions (alfa0, alfa1, ..., alfa(k - 1))
# If u(a) = u0, then u0 = (alfa0(p), alfa1(p), ..., alfa(k - 1)(p), p0, p1, ..., p(n-1))
# p - vector of last n - k components of vector u0
# Dimension of p is n - k
# Those expression is equal border conditions G(u(a)) = 0
# Dimension of alfa is k
# D - vector of functions from border conditions D(u(b)) = 0
# Dimension of D is n - k
# p0 - first iteration of parameters p
# x - argument of searched value u(x)
# h - step of iteration
# e - precision
def borderProblem_byCauchy(F, n, a, b, alfa, D, k, x, p0, h, w, eps):
    def stitching(p):
        return D(CauchyProblemParameters(F, n, a, alfa, k, p, b, h))

    p = secantSol(stitching, n - k, p0, w, eps)
    u0 = [0 for i in range(n)]

    for i in range(n):
        if i < k:
            u0[i] = alfa(p)[i]
        else:
            u0[i] = p[i - k]

    return rungeKutta4(n, F, a, u0, x, h)

# Example of using
#
# y''(x) + y'(x) exp(x+2) + y(x) sin(x) = tg(x+1)
# y(0) = 0
# y(1) = 0
#
# u(x) = (y(x), y'(x))
# f(x, u) = (u2, tan(x+1) - u2 exp(x + 2) - u1 sin(x))
# u'(x) = f(x, u(x))
#
# y(0) <=> u(0) = (alfa(p), p), alfa(p) = 0
# y(1) <=> D(u) = (u1 - 1)

def f(x, u):
    return [u[1], math.tan(x+1) - u[1] * math.exp(x + 2) - u[0] * math.sin(x)]

xa = 0
ya = 0

xb = 1
yb = 1

def alfa(p):
    return [ya]

def D(u):
    return [u[0] - yb]

def beta(p):
    return [yb]

def G(u):
    return [u[0] - ya]

N = 51

w = 0.0000001
eps = 0.0000001

h = (xb - xa) / (N - 1)
X = [xa + i * h for i in range(N)]

for i in range(N):
    print(borderProblem_byCauchy(f, 2, xa, xb, alfa, D, 1, X[i], [0], h, w, eps)[0])

print("\n")

for i in range(N):
    print(borderProblem_byCauchy(f, 2, xb, xa, beta, G, 1, X[i], [0], h, w, eps)[0])
