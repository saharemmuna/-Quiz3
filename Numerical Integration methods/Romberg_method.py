import numpy as np
from colors import bcolors
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from sympy.utilities.lambdify import lambdify


def romberg_integration(func, a, b, n):
    """
    Romberg Integration

    Parameters:
    func (function): The function to be integrated.
    a (float): The lower limit of integration.
    b (float): The upper limit of integration.
    n (int): The number of iterations (higher value leads to better accuracy).

    Returns:
    float: The approximate definite integral of the function over [a, b].
    """
    h = b - a
    R = np.zeros((n, n), dtype=float)

    R[0, 0] = 0.5 * h * (func(a) + func(b))

    for i in range(1, n):
        h /= 2
        sum_term = 0

        for k in range(1, 2 ** i, 2):
            sum_term += func(a + k * h)


        R[i, 0] = 0.5 * R[i - 1, 0] + h * sum_term

        for j in range(1, i + 1):
            R[i, j] = R[i, j - 1] + (R[i, j - 1] - R[i - 1, j - 1]) / ((4 ** j) - 1)
            print("i=", i - 1, "j=", j - 1, R[i - 1, j - 1])

    return R[n - 1, n - 1]

def calc_error(a, b, n):
    h = (b-a)/n
    x = sp.symbols('x')
    # Define the function f using SymPy
    f = x ** 3 + 2 * x ** 2 - x + 1
    print(f"f: {f}")

    # First derivative of f with respect to x
    df = sp.diff(f, x)
    print(f"df: {df}")

    # Second derivative of f with respect to x
    ddf = sp.diff(df, x)
    print(f"ddf: {ddf}")
    ddf = lambdify(x, ddf)

    x0 = 2 # point in which you want to find error

    E_x0 = (1/12) * (h**2) * (b-a) * ddf(x0)
    print(f"Error in {x0} is: " + str(E_x0))

def f(x):
    return 1/(2+x ** 4)


if __name__ == '__main__':

    a = 0
    b = 1
    n = 5
    integral = romberg_integration(f, a, b, n)

    print( f" Division into n={n} sections ")
    print(bcolors.OKBLUE, f"Approximate integral in range [{a},{b}] is {integral}", bcolors.ENDC)

