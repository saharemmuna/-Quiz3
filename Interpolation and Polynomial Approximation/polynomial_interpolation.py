from colors import bcolors
from matrix_utility import *
from jacobi_utilities import *

import numpy as np
def forward_substitution_to_diagonal(mat):
    N = len(mat)
    for k in range(N - 1, -1, -1):
        scalar = mat[k][k]
        for j in range(N + 1):
            mat[k][j] /= scalar

        for i in range(k - 1, -1, -1):
            scalar = mat[i][k]
            for j in range(N + 1):
                mat[i][j] -= mat[k][j] * scalar

def gaussianElimination(mat):
    N = len(mat)
    singular_flag = forward_substitution(mat)
    if singular_flag != -1:

        if mat[singular_flag][N]:
            return "Singular Matrix (Inconsistent System)"
        else:
            return "Singular Matrix (May have infinitely many solutions)"

# if matrix is non-singular:
    forward_substitution_to_diagonal(mat)
    #print(np.array(mat))
    # get solution to system using backward substitution
    return backward_substitution(mat)

def polynomialInterpolation(matrix,table_points, x):

    b = [[point[1]] for point in table_points]
    matrixNew = np.hstack((matrix, b))
    #print(bcolors.OKBLUE, "The matrix obtained from the points: ", bcolors.ENDC,'\n', np.array(matrix))
    #print(bcolors.OKBLUE, "\nb vector: ", bcolors.ENDC,'\n',np.array(b))
    matrixSol = gaussianElimination(matrixNew)
    if matrixSol is not None:
        #print(bcolors.OKBLUE, "\nResult Gauss: ", bcolors.ENDC, '\n', np.array(matrixSol))
        result = sum([matrixSol[i] * (x ** i) for i in range(len(matrixSol))])
        #print(bcolors.OKBLUE, "\nThe polynom:", bcolors.ENDC)
        #print('P(X) = '+'+'.join([ '('+str(matrixSol[i])+') * x^' + str(i) + ' ' for i in range(len(matrixSol))]))
        print(bcolors.OKGREEN, f"\nThe Result of P(X={x}) is:", bcolors.ENDC)
        print(result)
        return result
    return None

def Prerequisite(table_points):
    matrix = [[point[0] ** i for i in range(len(table_points))] for point in table_points]  # Makes the initial matrix
    if not np.linalg.det(matrix):
        print("Singular Matrix")
        return None
    return matrix

if __name__ == '__main__':
    """"
               Date: 8/4/24
               Group: Avishag Tamssut id-326275609
                       Merav Hashta id-214718405
                       Sahar Emmuna id-213431133
               Git: https://github.com/Avishagtams/Numerical-Analysis-Quiz2.git
               Name: Sahar Emmuna id-213431133
               """
    print(bcolors.OKBLUE, "----------------- Interpolation & Extrapolation Methods -----------------", bcolors.ENDC)
    table_points = [(0.2, 3.7241), (0.35, 3.9776), (0.45, 4.0625), (0.6, 2.9776), (0.75, 3.7241)]
    x = 0.4
    x2 = 0.65
    matrix = Prerequisite(table_points)
    if matrix is not None:
        #print(bcolors.OKBLUE, "Table Points: ", bcolors.ENDC, table_points)
        #print(bcolors.OKBLUE, "Finding an approximation to the point: ", bcolors.ENDC, x,'')
        if polynomialInterpolation(matrix, table_points, x) is None:
            print(" Singular Matrix")
        print(bcolors.OKBLUE, "---------------------------------------------------------------------------", bcolors.ENDC)
        if polynomialInterpolation(matrix, table_points, x2) is None:
            print(" Singular Matrix")
        print(bcolors.OKBLUE, "---------------------------------------------------------------------------",
              bcolors.ENDC)

'''
from colors import bcolors
from matrix_utility import *


def GaussJordanElimination(matrix, vector):
    """
    Function for solving a linear equation using gauss's elimination method
    :param matrix: Matrix nxn
    :param vector: Vector n
    :return: Solve Ax=b -> x=A(-1)b
    """
    # Pivoting process
    matrix, vector = RowXchange(matrix, vector)
    # Inverse matrix calculation
    invert = InverseMatrix(matrix, vector)
    return MulMatrixVector(invert, vector)




def UMatrix(matrix,vector):
    """
    :param matrix: Matrix nxn
    :return:Disassembly into a  U matrix
    """
    # result matrix initialized as singularity matrix
    U = MakeIMatrix(len(matrix), len(matrix))
    # loop for each row
    for i in range(len(matrix[0])):
        # pivoting process
        matrix, vector = RowXchageZero(matrix, vector)
        for j in range(i + 1, len(matrix)):
            elementary = MakeIMatrix(len(matrix[0]), len(matrix))
            # Finding the M(ij) to reset the organs under the pivot
            elementary[j][i] = -(matrix[j][i])/matrix[i][i]
            matrix = MultiplyMatrix(elementary, matrix)
    # U matrix is a doubling of elementary matrices that we used to reset organs under the pivot
    U = MultiplyMatrix(U, matrix)
    return U


def LMatrix(matrix, vector):
    """
       :param matrix: Matrix nxn
       :return:Disassembly into a  L matrix
       """
    # Initialize the result matrix
    L = MakeIMatrix(len(matrix), len(matrix))
    # loop for each row
    for i in range(len(matrix[0])):
        # pivoting process
        matrix, vector = RowXchageZero(matrix, vector)
        for j in range(i + 1, len(matrix)):
            elementary = MakeIMatrix(len(matrix[0]), len(matrix))
            # Finding the M(ij) to reset the organs under the pivot
            elementary[j][i] = -(matrix[j][i])/matrix[i][i]
            # L matrix is a doubling of inverse elementary matrices
            L[j][i] = (matrix[j][i]) / matrix[i][i]
            matrix = MultiplyMatrix(elementary, matrix)

    return L


def SolveLU(matrix, vector):
    """
    Function for deconstructing a linear equation by ungrouping LU
    :param matrix: Matrix nxn
    :param vector: Vector n
    :return: Solve Ax=b -> x=U(-1)L(-1)b
    """
    matrixU = UMatrix(matrix)
    matrixL = LMatrix(matrix)
    return MultiplyMatrix(InverseMatrix(matrixU), MultiplyMatrix(InverseMatrix(matrixL), vector))


def solveMatrix(matrixA,vectorb):
    detA = Determinant(matrixA, 1)
    print(bcolors.YELLOW, "\nDET(A) = ", detA)

    if detA != 0:
        print("CondA = ", Cond(matrixA, InverseMatrix(matrixA, vectorb)), bcolors.ENDC)
        print(bcolors.OKBLUE, "\nnon-Singular Matrix - Perform GaussJordanElimination",bcolors.ENDC)
        result = GaussJordanElimination(matrixA, vectorb)
        print(np.array(result))
        return result
    else:
        print("Singular Matrix - Perform LU Decomposition\n")
        print("Matrix U: \n")
        print(np.array(UMatrix(matrixA, vectorb)))
        print("\nMatrix L: \n")
        print(np.array(LMatrix(matrixA, vectorb)))
        print("\nMatrix A=LU: \n")
        result = MultiplyMatrix(LMatrix(matrixA, vectorb), UMatrix(matrixA, vectorb))
        print(np.array(result))
        return result


def polynomialInterpolation(table_points, x):
    matrix = [[point[0] ** i for i in range(len(table_points))] for point in table_points] # Makes the initial matrix

    b = [[point[1]] for point in table_points]

    print(bcolors.OKBLUE, "The matrix obtained from the points: ", bcolors.ENDC,'\n', np.array(matrix))
    print(bcolors.OKBLUE, "\nb vector: ", bcolors.ENDC,'\n',np.array(b))
    matrixSol = solveMatrix(matrix, b)

    result = sum([matrixSol[i][0] * (x ** i) for i in range(len(matrixSol))])
    print(bcolors.OKBLUE, "\nThe polynom:", bcolors.ENDC)
    print('P(X) = '+'+'.join([ '('+str(matrixSol[i][0])+') * x^' + str(i) + ' ' for i in range(len(matrixSol))])  )
    print(bcolors.OKGREEN, f"\nThe Result of P(X={x}) is:", bcolors.ENDC)
    print(result)
    return result


if __name__ == '__main__':

    table_points = [(1, 3), (2, 4), (3, -1)]
    x = 1.5
    print(bcolors.OKBLUE, "----------------- Interpolation & Extrapolation Methods -----------------\n", bcolors.ENDC)
    print(bcolors.OKBLUE, "Table Points: ", bcolors.ENDC, table_points)
    print(bcolors.OKBLUE, "Finding an approximation to the point: ", bcolors.ENDC, x,'\n')
    polynomialInterpolation(table_points, x)
    print(bcolors.OKBLUE, "\n---------------------------------------------------------------------------\n", bcolors.ENDC)
'''
