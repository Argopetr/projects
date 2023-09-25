import math
# a(x)y''(x) + b(x)y'(x) + c(x)y(x) = f(x),  x c [xa, xb]
# y(xa) = ya
# y(xb) = yb
# N - number of points
# w - parametre of linear iteration
# eps - precision
def boardProblem2_Differents(a, b, c, f, xa, xb, ya, yb, N):
    h = (xb - xa) / (N - 1)
    
    A = [a(xa + i * h) / (h * h) - b(xa + i * h) / (2 * h) for i in range(N)]
    B = [c(xa + i * h) - 2 * a(xa + i * h) / (h * h) for i in range(N)]
    C = [a(xa + i * h) / (h * h) + b(xa + i * h) / (2 * h) for i in range(N)]
    F = [f(xa + i * h) for i in range(N)]
    alfa = [0 for i in range(N)]
    beta = [0 for i in range(N)]
    X = [xa + i * h for i in range(N)]
    Y = [0 for i in range(N)]

    F[1] -= A[1] * ya
    A[1] = 0

    F[N - 2] -= C[N - 2] * yb
    C[N - 2] = 0

    for i in range(1, N - 1):
        alfa[i + 1] = - C[i] / (A[i] * alfa[i] + B[i])
        beta[i + 1] = (F[i] - A[i] * beta[i]) / (A[i] * alfa[i] + B[i])

    Y[0] = ya
    Y[N - 2] = (F[N - 2] - A[N - 2] * beta[N - 2]) / (A[N - 2] * alfa[N - 2] + B[N - 2])
    Y[N - 1] = yb

    for i in range(N - 3):
        Y[N - 3 - i] = alfa[N - 2 - i] * Y[N - 2 - i] + beta[N - 2 - i]

    return [X, Y]

def a(x):
    return 1

def b(x):
    return math.exp(x + 2)

def c(x):
    return math.sin(x)

def f(x):
    return math.tan(x + 1)

xa = 0
ya = 0

xb = 1
yb = 1

N = 51

h = (xb - xa) / (N - 1)
X = [xa + i * h for i in range(N)]

A = boardProblem2_Differents(a, b, c, f, xa, xb, ya, yb, N)
B = boardProblem2_Differents(a, b, c, f, xb, xa, yb, ya, N)

for i in range(N):
    print(A[0][i], A[1][i])

print()

for i in range(N):
    print(B[0][N - 1 - i], B[1][N - 1 - i])

