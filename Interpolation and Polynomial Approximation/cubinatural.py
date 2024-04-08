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


def spline(f, x0, Fx_0,Fx_n):
    n = len(f)

    matrix = [[0] * (n+1) for _ in range(n)]

    h0 = f[1][0]-f[0][0]
    h_n_1 = f[n-1][0]-f[n-2][0]

    matrix[0][0] = h0/3
    matrix[0][1] = h0/6
    matrix[n-1][n-1] = h_n_1/3
    matrix[n-1][n-2] = h_n_1/6

    matrix[0][n] = (f[1][1]-f[0][1])/h0 - Fx_0
    matrix[n - 1][n] = Fx_n - (f[n-1][1]-f[n-2][1])/h_n_1



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

    print_xo(f,S_list, x0)







def calc_im_spline():
    f = [(1, 1), (2, 2), (3, 1), (4, 1.5), (5, 1)]
    h = [(-1, 1), (0, 0), (1, 1)]
    x0 = 0.5
    spline(h, x0, -1, -1)



if __name__ == '__main__':
    calc_im_spline()
