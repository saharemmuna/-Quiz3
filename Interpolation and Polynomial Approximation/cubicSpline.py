import numpy as np
import sympy as sp
from jacobi_utilities import *
from colors import bcolors
x = sp.symbols('x')


def print_xo(f, S_list, x0):
    n = len(f)
    if (x0 < f[0][0]):
        print(
            "\nx0 smaller than f(x" + str(0) + ") = " + str(f[0][0]) + " so:")
        print("s" + str(0) + "(" + str(x0) + ") = " + str(float(S_list[0].subs(x, x0))))
    else:
        if (x0 > f[n - 1][0]):
            print(
                "\nx0 bigger than f(x" + str(n) + ") = " + str(f[n - 1][0]) + " so:")
            print("s" + str(n - 1) + "(" + str(x0) + ") = " + str(float(S_list[n - 2].subs(x, x0))))

        else:
            for i in range(n - 1):
                if (x0 > f[i][0] and x0 < f[i + 1][0]):
                    print(
                        "\nx0 between f(x" + str(i + 1) + ") = " + str(f[i][0]) + " and f(x" + str(
                            i + 2) + ") = " + str(
                            f[i + 1][0]) + " so:")
                    print("s" + str(i + 1) + "(" + str(x0) + ") = " + str(float(S_list[i].subs(x, x0))))


def spline(f, x0):
    n = len(f)

    matrix = [[0] * (n+1) for _ in range(n)]

    matrix[0][0] = 1
    matrix[n-1][n-1] = 1



    for i in range(1, n-1):
        h0 = f[i][0]-f[i-1][0]
        h1 = f[i+1][0]-f[i][0]

        matrix[i][n] = (f[i+1][1]-f[i][1])/h1 - (f[i][1]-f[i-1][1])/h0

        for j in range(0, n):
            if i == j:
                matrix[i][j] = (h0+h1)/3
                matrix[i][j-1] = h0/6
                matrix[i][j+1] = h1/6

    print(np.array(matrix))

    solution = gaussianElimination(matrix)

    if isinstance(solution, str):
        print(solution)
    else:
        print(bcolors.OKBLUE, "\nSolution for the system:")
        for x1 in solution:
            print("{:.6f}".format(x1))



    S_list = [None] * n

    for i in range(n-1):
        h_i = f[i+1][0] - f[i][0]
        S_list[i] = ((f[i+1][1] * (x-f[i][0]))/h_i - (f[i][1] * (x-f[i+1][0]))/h_i
                               + solution[i+1]/6 * (((x-f[i][0]) ** 3)/h_i - h_i * (x-f[i][0])) -
                               solution[i]/6 * (((x-f[i+1][0]) ** 3)/h_i - h_i * (x-f[i+1][0])))

        print("s" + str(i) + "(x) = " + str(S_list[i]))

    print_xo(f,S_list,x0)







def calc_im_spline():
    f = [(1, 1), (2, 2), (3, 1), (4, 1.5), (5, 1)]
    x0 = 4.5
    spline(f, 4)



if __name__ == '__main__':
    calc_im_spline()




"""
import jacobi_utilities
from sympy import *
import numpy as np

x = Symbol('x')


def natural_cubic_spline(f, x0):
    h = list()
    for i in range(len(f) - 1):
        h.append(f[i + 1][0] - f[i][0])

    g = list()
    g.append(0)  # g0
    for i in range(1, len(f) - 1):
        g.append(h[i] / (h[i] + h[i - 1]))
    g.append(0)  # gn

    m = list()
    m.append(0)
    for i in range(1, len(f)):
        m.append(1 - g[i])


    d = list()
    d.append(0)  # d0=0
    for i in range(1, len(f) - 1):
        d.append((6 / (h[i - 1] + h[i])) * (((f[i + 1][1] - f[i][1]) / h[i]) - ((f[i][1] - f[i - 1][1]) / h[i - 1])))
    d.append(0)  # dn

    # building the matrix
    mat = list()

    # first row
    mat.append(list())
    mat[0].append(2)
    for j in range(len(f) - 1):
        mat[0].append(0)

    for i in range(1, len(f) - 1):
        mat.append(list())
        for j in range(len(f)):
            if j == i - 1:  # put miu
                mat[i].append(m[i])
            elif j == i:
                mat[i].append(2)
            elif j == i + 1:  # put lambda
                mat[i].append(g[i])
            else:
                mat[i].append(0)

    # last row
    mat.append(list())
    for j in range(len(f) - 1):
        mat[len(f) - 1].append(0)
    mat[len(f) - 1].append(2)

    print("matrix: " + str(mat))
    print("vector b: " + str(d))

    # get m vector
    print("\nJacobi middle results: ")
    M = (jacobi_utilities.Jacobi(mat, d))
    print("\nvector M: " + str(list(map(float, M))))

    # find S:
    for loc in range(1, len(f)):
        s = (((f[loc][0] - x) ** 3) * M[loc - 1] + ((x - f[loc - 1][0]) ** 3) * M[loc]) / (6 * h[loc - 1])
        s += (((f[loc][0] - x) * f[loc - 1][1]) + ((x - f[loc - 1][0]) * f[loc][1])) / h[loc - 1]
        s -= (((f[loc][0] - x) * M[loc - 1] + (x - f[loc - 1][0]) * M[loc]) * h[loc - 1]) / 6
        print("s" + str(loc - 1) + "(x) = " + str(s))

    # find the location of x0:
    loc = 0
    for i in range(1, len(f)):
        if x0 < f[i][0] and x0 > f[i - 1][0]:
            loc = i
            break

    if loc == 0:
        print("no range found for x0")
        return

    s = (((f[loc][0] - x) ** 3) * M[loc - 1] + ((x - f[loc - 1][0]) ** 3) * M[loc]) / (6 * h[loc - 1])
    s += (((f[loc][0] - x) * f[loc - 1][1]) + ((x - f[loc - 1][0]) * f[loc][1])) / h[loc - 1]
    s -= (((f[loc][0] - x) * M[loc - 1] + (x - f[loc - 1][0]) * M[loc]) * h[loc - 1]) / 6

    print("\nx0 between f(x" + str(loc - 1) + ") = " + str(f[loc - 1][0]) + " and f(x" + str(loc) + ") = " + str(
        f[loc][0]) + " so:")
    print("s" + str(loc - 1) + "(" + str(x0) + ") = " + str(float(s.subs(x, x0))))

def natural_cubics_spline(f, x0, df=None):
    lengh = len(f)
    h = list()
    for i in range(len(f) - 1):
        h.append(f[i + 1][0] - f[i][0])
    A = np.zeros((len(f), len(f)))
    b=np.zeros(len(f))
    for i in range (len(f)):
        for j in range (len(f)):
            if i==0 :
                if df is not None:
                    A[0][0] = 1/3 * h[0]
                    A[0][1] = 1/6 * h[0]
                    b[0] = (f[1][1]-f[0][1]) / h[0] - df[0]
                else:
                    A[0][0] = 1
                    b[0] =0


            if i == len(f)-1 :
                if df is not None:
                    A[i][i-1] = 1/6 * h[i-1]
                    A[i][i] = 1/3 * h[i-1]
                    b[i] = df[1] - (f[i][1]-f[i-1][1]) / h[i-1]













if __name__ == '__main__':
    f = [(1, 1), (2, 0), (5, 2)]

    x0 = 3


    print("func: " + str(f))
    print("x0 = " + str(x0) + "\n")
    natural_cubic_spline(f, x0)

"""
"""
from jacobi_utilities import *
from sympy import *
from matrix_utility import *
import numpy as np

x = Symbol('x')




def CubicSpline(tableValue, X, Ftag = None):
    size = len(tableValue)
    matrix = np.zeros((size, size))
    b = np.zeros(size)
    h = [tableValue[i+1][0] - tableValue[i][0] for i in range(size - 1)]
    if Ftag == None:
        matrix[0][0] = 1
        b[0] = 0
        matrix[-1][-1] = 1
        b[-1] = 0
    else:
        matrix[0][0] = 1/3 * h[0]
        matrix[0][1] = 1/6 * h[0]
        b[0] = (tableValue[1][1] - tableValue[0][1]) - Ftag[0]
        matrix[-1][-1] = 1/6 * h[-1]
        matrix[-1][-2] = 1/3 * h[-1]
        b[-1] = Ftag[1] - ((tableValue[-1][1] - tableValue[-2][1])/ h[-1])
    for i in range(1, size - 1):
        for j in range(i - 1, min(i + 2, size)):
            if i != j:
                matrix[i][j] = 1/6 * h[i - 1]
            else:
                matrix[i][j] = 1/3 * (h[i - 1] + h[i])
        b[i] = ((tableValue[i+1][1] - tableValue[i][1]) / h[i]) - ((tableValue[i][1] - tableValue[i-1][1])/h[i-1])
    matrixNew = np.hstack((matrix, b.reshape(-1, 1)))
    matrixSol = gaussianElimination(matrixNew)
    sum = 0
    S_Func = 0
    x = Symbol('x')
    print(np.array(matrixSol))
    for i in range(size - 1):
        #print(tableValue[i+1][1])
        #print(matrixSol[i + 1])
        sum += ((tableValue[i+1][1] * (x - tableValue[i][0]))/h[i])
        sum -= (tableValue[i][1] * (x - tableValue[i+1][0])/h[i])
        sum += ((matrixSol[i+1]/6) * ((((x - tableValue[i][0])**3)/h[i]) - h[i] * (x - tableValue[i][0])))
        sum -= (matrixSol[i]/6) * ((((x - tableValue[i+1][0])**3)/h[i]) - (h[i] * (x - tableValue[i+1][0])))
        print("s" + str(i) + "(x) = " + str(sum))
        S_Func += sum
        sum = 0
    print("S=", S_Func)
    s = lambdify(x, S_Func)
    print("------The sol------")
    print(s(3))





if __name__ == '__main__':
    f = [(1, 1), (2, 2), (3, 1), (4, 1.5), (5, 1)]
    x0 = 2
    print("func: " + str(f))
    print("x0 = " + str(x0) + "\n")
    CubicSpline(f, x0)
"""
