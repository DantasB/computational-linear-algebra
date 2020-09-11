import copy
import math

def display_matrix(matrix):
    print('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in matrix]))

    print('\n')


def matrix_determinant(matrix):
    number_of_rows    = len(matrix)
    number_of_columns = len(matrix[0])
    result  = 0

    if(number_of_rows != number_of_columns):
        return "Error"

    if(number_of_rows == 1):
        return matrix[0][0]

    for k in range(number_of_rows):
        result += matrix[k][0]*((-1)**k)*matrix_determinant(get_auxiliar_matrix(matrix, k))

    return result


def get_auxiliar_matrix(matrix, index):
    secondary = copy.deepcopy(matrix)

    for row in range(len(secondary)):
        secondary[row] = secondary[row][1:]

    return secondary[:index] + secondary[index+1:]


def positive_definite(matrix):
    if(not check_simetry(matrix)):
        return False
    if(not sylvesters_criterion(matrix)):
        return False

    return True


def check_simetry(matrix):
    number_of_rows    = len(matrix)
    number_of_columns = len(matrix[0])

    if(number_of_rows != number_of_columns):
        return -1

    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if(matrix[i][j]!= matrix[j][i]):
                return False

    return True


def sylvesters_criterion(matrix):
    for i in range(len(matrix)):
        secondary = get_main_minor(matrix, i)
        if(matrix_determinant(secondary) <= 0):
            return False

    return True


def get_main_minor(matrix, index):
    result = [[matrix[row][column] for row in range(index + 1)] for column in range(index + 1)]
    return result


def backward_substitution(matrix_u, matrix_y):
    number_of_rows = len(matrix_u)
    matrix_x = [0 for i in range(number_of_rows)]

    matrix_x[number_of_rows-1] = matrix_y[number_of_rows-1]/matrix_u[number_of_rows-1][number_of_rows-1]

    for i in range(number_of_rows-2, -1, -1):
        summation = matrix_y[i]
        for j in range(i+1, number_of_rows):
            summation -= matrix_u[i][j]*matrix_x[j]

        matrix_x[i] = summation/float(matrix_u[i][i])
    return matrix_x


def forward_substitution(matrix_l, matrix_b, control=False):
    number_of_rows = len(matrix_l)
    matrix_y = [0 for i in range(number_of_rows)]

    if(not control):
        matrix_y[0] = matrix_b[0]
    else:
        matrix_y[0] = matrix_b[0]/matrix_l[0][0]

    for i in range(1, number_of_rows):
        summation = matrix_b[i]
        for j in range(i):
            summation -= matrix_l[i][j]*matrix_y[j]

        if(not control):
            matrix_y[i] = summation
        else:
            matrix_y[i] = summation/matrix_l[i][i]

    return matrix_y


def converge(matrix):
    for i in range(len(matrix)):
        lines_summation   = 0
        columns_summation = 0;
        for j in range(len(matrix)):
            if (i != j):
                lines_summation+= math.fabs(matrix[i][j])
                columns_summation+=math.fabs(matrix[j][i])

        if(matrix[i][i] < lines_summation or matrix[i][i] < columns_summation):
            return False
    return True


def sum_of_vector_multiplication(matrix_a, matrix_b):
    summation = 0

    for i in range(len(matrix_a)):
        for j in range(len(matrix_b)):
            if (i==j):
                summation += matrix_a[i]*matrix_b[i]

    return summation


def get_transposed_matrix(matrix):
    answer = [[0.0]*len(matrix) for i in range(len(matrix))]

    for i in range(len(matrix)):
        for j in range(len(matrix)):
            answer[j][i] = matrix[i][j]

    return answer
