
import math
import numpy as np
from colors import bcolors

"""
Receives 3 parameters:
    1.  a - start value.
    2.  b - end  value. 
    3.  err - value of tolerable error

Returns variables:
    1.  S - The minimum number of iterations required to reach the desired accuracy
"""


def max_steps(a, b, err):
    s = int(np.floor(- np.log2(err / (b - a)) / np.log2(2) - 1))
    return s


"""
Performs Iterative methods for Nonlinear Systems of Equations to determine the roots of the given function f
Receives 4 parameters:
    1. f - continuous function on the interval [a, b], where f (a) and f (b) have opposite signs.
    2. a - start value.
    3. b - end  value. 
    4. tol - the tolerable error , the default value will set as 1e-16

Returns variables:
    1.  c - The approximate root of the function f
"""


def bisection_method(f, a, b, tol=1e-6):
    # if np.sign(a) == np.sign(b):
    #     raise Exception("The scalars a and b do not bound a root")
    c, k = 0, 0
    steps = max_steps(a, b, tol)  # calculate the max steps possible

    print("{:<10} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}".format("Iteration", "a", "b", "f(a)", "f(b)", "c", "f(c)"))

    # while the diff af a&b is not smaller than tol, and k is not greater than the max possible steps
    while abs(b - a) > tol and k <= steps:
        c = (a + b) / 2  # Calculation of the middle value

        if f(c) == 0:
            return c  # Procedure completed successfully

        if f(c) * f(a) < 0:  # if sign changed between steps
            b = c  # move forward
        else:
            a = c  # move backward

        print("{:<10} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f}".format(k, a, b, f(a), f(b), c, f(c)))
        k += 1

    return c  # return the current root


def find_all_roots(f, a, b, tol=1e-6):
    roots = []
    x = np.linspace(a, b, 1000)  # Divide the interval into smaller sub-intervals

    for i in range(len(x) - 1):
        if np.sign(f(x[i])) != np.sign(f(x[i + 1])):
            root = bisection_method(f, x[i], x[i + 1], tol)
            roots.append(root)

    return roots


if __name__ == '__main__':
    """"
               Date: 18/3/24
               Group: Avishag Tamssut id-326275609
                       Merav Hashta id-214718405
                       Sahar Emmuna id-213431133
               Git: https://github.com/Avishagtams/Numerical-Analysis-Quiz2.git
               Name: Avishag Tamssut 326275609
               """
    f = lambda x: (5*x**3 + 1*x**2 + 2) / (2*x - 5)
    a = -3
    b = 0

    roots = find_all_roots(f, a, b)
    print(bcolors.OKBLUE, f"\nThe equation f(x) has approximate roots at {roots}", bcolors.ENDC)




'''
import math
import numpy as np
from colors import bcolors

"""
Receives 3 parameters:
    1.  a - start value.
    2.  b - end  value. 
    3.  err - value of tolerable error

Returns variables:
    1.  S - The minimum number of iterations required to reach the desired accuracy
"""
def max_steps(a, b, err):
    s = int(np.floor(- np.log2(err / (b - a)) / np.log2(2) - 1))
    return s

"""
Performs Iterative methods for Nonlinear Systems of Equations to determine the roots of the given function f
Receives 4 parameters:
    1. f - continuous function on the interval [a, b], where f (a) and f (b) have opposite signs.
    2. a - start value.
    3. b - end  value. 
    4. tol - the tolerable error , the default value will set as 1e-16

Returns variables:
    1.  c - The approximate root of the function f
"""
def bisection_method(f, a, b, tol=1e-6):
    if np.sign(f(a)) == np.sign(f(b)):
        raise Exception("The scalars a and b do not bound a root")
    c, k = 0, 0
    steps = max_steps(a, b, tol)  # calculate the max steps possible

    print("{:<10} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}".format("Iteration", "a", "b", "f(a)", "f(b)", "c", "f(c)"))

    # while the diff af a&b is not smaller than tol, and k is not greater than the max possible steps
    while abs(b - a) > tol and k < steps:
        c = a + (b - a) / 2  # Calculation of the middle value

        if f(c) == 0 :
            return c  # Procedure completed successfully

        if f(c) * f(a) < 0:  # if sign changed between steps
            b = c  # move forward
        else:
            a = c  # move backward

        print("{:<10} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f}".format(k, a, b, f(a), f(b), c, f(c)))
        k += 1

    return c  # return the current root

def iterative_method(f, x0, tol=1e-6, max_iter=100):
    """
    פונקציה זו מממשת שיטת איטרציות פשוטה לפתרון משוואות באמצעות הנוסחה Xr+1=f(Xr)+Xr.

    :param f: הפונקציה f(x).
    :param x0: ערך ראשוני כלשהו ל־x.
    :param tol: דיוק מבוקש לתוצאה.
    :param max_iter: מספר מקסימלי של איטרציות.
    :return: הערך הקרוב ביותר לנקודת האפס.
    """
    x = x0
    for k in range(max_iter):
        x_new = f(x) + x
        if abs(x_new - x) < tol:
            #print(f"Iteration {k}: x = {x}, f(x) = {f(x)}, x_new = {x_new}")
            return x_new
        #print(f"Iteration {k}: x = {x}, f(x) = {f(x)}, x_new = {x_new}")
        x = x_new
    return None  # לא מצליח למצוא אפס תוך מספר האיטרציות המקסימל


if __name__ == '__main__':
    f = lambda x: math.cos(x)  # פונקציית הדוגמה: f(x) = cos(x)
    roots = bisection_method(f, 1, 3)
    print(bcolors.OKBLUE, f"\nThe equation f(x) has an approximate root at x = {roots}",bcolors.ENDC,)


f = lambda x: math.cos(x)  # פונקציית הדוגמה: f(x) = cos(x)
root_it = iterative_method(f, 0)  # מחפשים את האפס כשהתחלת האיטרציה היא ב־0.5
print(bcolors.OKBLUE, f"\nIterative method equation f(x) has an approximate root at x = {root_it}", bcolors.ENDC, )

#print("f(root):", f(root))

'''