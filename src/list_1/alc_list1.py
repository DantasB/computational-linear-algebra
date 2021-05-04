import copy
import math
from src.utils.matrix_utils import vector_multiplication, converge, matrix_determinant, positive_definite, forward_substitution, backward_substitution, get_transposed_matrix


def lu_decomposition(matrix):
    if(matrix_determinant(matrix) == 0):
        return "Error"

    number_of_rows = len(matrix)
    number_of_columns = len(matrix[0])

    if(number_of_rows != number_of_columns):
        return "Error"

    result = copy.deepcopy(matrix)

    for k in range(number_of_rows):
        for i in range(k+1, number_of_rows):
            result[i][k] = float(result[i][k]/result[k][k])

        for j in range(k+1, number_of_columns):
            for i in range(k+1, number_of_columns):
                result[i][j] = float(result[i][j]-result[i][k]*result[k][j])

    return result


def cholesky_decomposition(matrix):
    number_of_rows = len(matrix)
    number_of_columns = len(matrix[0])

    if(number_of_rows != number_of_columns):
        return "Error"

    # Check if the matrix is simetric, if not, return Error
    if(not positive_definite(matrix)):
        return "Error"

    l = [[0.0] * len(matrix) for _ in range(len(matrix))]

    for i in range(len(matrix)):
        for j in range(i + 1):

            if(i == j):
                summation = sum(l[i][k]**2 for k in range(i))
                l[i][i] = (matrix[i][i]-summation)**0.5
                continue

            summation = sum(l[i][k]*l[j][k] for k in range(i))
            l[i][j] = (1.0/l[j][j])*(matrix[i][j]-summation)

    return l


def solve_matrix(matrix_a, matrix_b, use_cholesky):
    matrix_lu = copy.deepcopy(matrix_a)
    if(not use_cholesky):
        matrix_lu = lu_decomposition(matrix_lu)
        if(matrix_lu == "Error"):
            return "Error"

        matrix_y = forward_substitution(matrix_lu, matrix_b)
        return backward_substitution(matrix_lu, matrix_y)

    matrix_lu = cholesky_decomposition(matrix_lu)
    if(matrix_lu == "Error"):
        return "Error"

    matrix_y = forward_substitution(matrix_lu, matrix_b, True)
    return backward_substitution(get_transposed_matrix(matrix_lu), matrix_y)


def iterative_jacobi(matrix_a, matrix_b):
    if (not converge(matrix_a)):
        return "Error"

    n = len(matrix_a)

    solution_zero = [0.0 for i in range(n)]
    next_solution = [0.0 for i in range(n)]

    tol = 10**(-5)
    residue = 1
    step = 0

    while (residue > tol):

        numerator = 0
        denominator = 0

        for j in range(n):
            next_solution[j] = matrix_b[j]

            for k in range(n):
                if (j != k):
                    next_solution[j] += (-1)*(matrix_a[j]
                                              [k] * solution_zero[k])

            next_solution[j] /= matrix_a[j][j]

        for z in range(n):
            numerator += (next_solution[z]-solution_zero[z])**2
            denominator += next_solution[z]**2

        residue = float(numerator**0.5)/(denominator**0.5)

        for i in range(len(next_solution)):
            solution_zero[i] = next_solution[i]
        step += 1

    print("x1: ", next_solution)
    print("residue: ", residue)
    print("iteration_number: ", step)


def gauss_seidel(matrix_a, matrix_b):
    if (not converge(matrix_a)):
        return "Error"

    n = len(matrix_a)

    solution_zero = [1.0 for i in range(n)]
    next_solution = [0.0 for i in range(n)]

    tol = 10**(-5)
    residue = 1
    step = 0

    while (residue > tol):
        numerator = 0
        denominator = 0
        second_summation = 0
        second_summation = 0

        for j in range(n):
            first_summation = vector_multiplication(
                matrix_a[j][:j], next_solution[:j])
            second_summation = vector_multiplication(
                matrix_a[j][j+1:], solution_zero[j+1:])
            next_solution[j] = (matrix_b[j] - first_summation -
                                second_summation)/matrix_a[j][j]

        for z in range(n):
            numerator += (next_solution[z] - solution_zero[z])**2
            denominator += next_solution[z]**2

        residue = float(numerator**0.5)/(denominator**0.5)

        for i in range(len(next_solution)):
            solution_zero[i] = next_solution[i]

        step += 1

    print("x1: ", next_solution)
    print("residue: ", residue)
    print("iteration_number: ", step)
