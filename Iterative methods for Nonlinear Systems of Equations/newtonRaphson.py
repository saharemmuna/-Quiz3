from colors import bcolors
import sympy as sp
from sympy.utilities.lambdify import lambdify

def newton_raphson(f, df, p0, TOL, N=50):
    print("{:<10} {:<15} {:<15} ".format("Iteration", "po", "p1"))
    for i in range(N):
        if df(p0) == 0:
            print( "Derivative is zero at p0, method cannot continue.")
            return None

        p = p0 - f(p0) / df(p0)

        if abs(p - p0) < TOL:
            return p  # Procedure completed successfully
        print("{:<10} {:<15.9f} {:<15.9f} ".format(i, p0, p))
        p0 = p
    return p


if __name__ == '__main__':
    f = lambda x: x**3 - 3*x**2
    df = lambda x: 3*x**2 - 6*x
    p0 = -3
    TOL = 1e-6
    N = 100
    f1 = lambda x: (5*x**3 + 1*x**2 + 2) / (2*x - 5)
    df1 = lambda x: (20 * x**3 - 73 * x**2 - 10 * x -4) / ((2 * x - 5)**2)


    roots = newton_raphson(f1, df1,p0,TOL,N)
    if roots is not None:
        print(bcolors.OKBLUE, "\nThe equation f(x) has an approximate root at x = {:<15.9f} ".format(roots),  bcolors.ENDC, )
    else:
        print("Unable to find root with the given parameters.")


'''
    x = sp.symbols('x')
    my_f = (5 * x ** 3 + 1 * x ** 2 + 2) / (2 * x - 5)
    print("my_func: ", my_f)  # my_func:  x**3 + 2*x + 5
    my_f1 = sp.diff(my_f, x)  # Derivation of my_f by x
    print("f' : ", my_f1)  # print my_f1
'''

