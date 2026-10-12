"""Small dense linear algebra for symmetric positive-definite matrices (pure Python).

Matrices are lists of row lists; vectors are lists.
"""

import math


class SingularMatrixError(ArithmeticError):
    pass


def cholesky(a, tol=1e-12):
    """Lower-triangular L with a = L L^T. Raises SingularMatrixError if a is
    not (numerically) positive definite."""
    p = len(a)
    scale = max((a[i][i] for i in range(p)), default=0.0)
    if scale <= 0:
        raise SingularMatrixError("matrix is not positive definite")
    L = [[0.0] * p for _ in range(p)]
    for j in range(p):
        d = a[j][j] - sum(L[j][k] ** 2 for k in range(j))
        if d <= tol * scale:
            raise SingularMatrixError("matrix is singular or not positive definite")
        L[j][j] = math.sqrt(d)
        for i in range(j + 1, p):
            L[i][j] = (a[i][j] - sum(L[i][k] * L[j][k] for k in range(j))) / L[j][j]
    return L


def forward_solve(L, b):
    """Solve L y = b for lower-triangular L."""
    y = []
    for i, row in enumerate(L):
        y.append((b[i] - sum(row[k] * y[k] for k in range(i))) / row[i])
    return y


def backward_solve(L, y):
    """Solve L^T x = y for lower-triangular L."""
    p = len(L)
    x = [0.0] * p
    for i in reversed(range(p)):
        x[i] = (y[i] - sum(L[k][i] * x[k] for k in range(i + 1, p))) / L[i][i]
    return x


def cho_solve(L, b):
    """Solve a x = b given the Cholesky factor L of a."""
    return backward_solve(L, forward_solve(L, b))


def inverse_spd(a):
    L = cholesky(a)
    p = len(a)
    cols = [cho_solve(L, [1.0 if i == j else 0.0 for i in range(p)]) for j in range(p)]
    return [[cols[j][i] for j in range(p)] for i in range(p)]


def quad_form_inv(a, v):
    """v^T a^-1 v for symmetric positive-definite a."""
    y = forward_solve(cholesky(a), v)
    return sum(t * t for t in y)


def submatrix(a, rows, cols):
    return [[a[i][j] for j in cols] for i in rows]


def matvec(a, v):
    return [sum(x * y for x, y in zip(row, v)) for row in a]
