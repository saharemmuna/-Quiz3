import math
from colors import bcolors
import sympy as sp
from sympy.utilities.lambdify import lambdify

def trapezoidal_rule(f, a, b, n):

    h = (b - a) / n
    T = f(a) + f(b)
    integral = 0.5 * T  # Initialize with endpoints

    for i in range(1, n):
        x_i = a + i * h
        integral += f(x_i)

    integral *= h

    return integral

def error_t(f_e, a, b, error,n):
    h = (b - a) / n
    x = sp.symbols('x')
    #print("my_func: ", f_e)  # my_func:  x**3 + 2*x + 5
    df = sp.diff(f_e, x)  # Derivation of my_f by x
    #print("f' : ", df)  # print my_f1
    ddf = sp.diff(df, x)
    #print("f'' : ", ddf)  # print my_f1
    f_tag_tag = lambdify(x, ddf)
    er=f_tag_tag(error)
    y = (1/12)*(b-a)*h**2
    finallError= y*er
    print("-----------------------------------------------------")
    print(bcolors.OKBLUE,"The Error: ", finallError,bcolors.ENDC)
    print("-----------------------------------------------------")







if __name__ == '__main__':
    f = lambda x: math.e ** (x ** 2 )
    result = trapezoidal_rule(f, 0, 1, 2)
    print("-----------------------------------------------------")
    print(bcolors.OKBLUE,"Approximate integral:", result, bcolors.ENDC)
    x = sp.symbols('x')
    f_e = math.e ** (x ** 2)
    error_t(f_e, 0, 1, 1, 2)

