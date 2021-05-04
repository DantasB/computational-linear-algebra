import copy
import math
import random


def initial_guess(n):
    result = [random.uniform(-1, 1) for _ in range(n)]
    print("Initial guess: " + str(result))
    return result


def step_forward_derivate(function, x, delta):
    numerator = function(x + delta) - function(x)
    denominator = delta
    result = numerator/denominator

    return result


def step_backward_derivate(function, x, delta):
    numerator = function(x) - function(x - delta)
    denominator = delta
    result = numerator/denominator

    return result


def central_derivate(function, x, delta):
    numerator = function(x + delta) - function(x - delta)
    denominator = 2*delta
    result = numerator/denominator

    return result


def get_jacobian_matrix(functions_list, first_solution):
    first_dimension = len(functions_list)
    second_dimension = len(first_solution)

    result = [[0 for _ in range(second_dimension)]
              for _ in range(first_dimension)]

    for i in range(first_dimension):
        for j in range(second_dimension):
            result[i][j] = partial_derivate(
                functions_list[i], first_solution, j)
    return result


def partial_derivate(multiple_variables_function, first_solution, derivate_index):
    delta = 10**(-10)
    aux = multiple_variables_function(first_solution)

    new_solution = first_solution[:]
    new_solution[derivate_index] += delta

    numerator = multiple_variables_function(new_solution) - aux
    denominator = delta
    result = numerator/denominator

    return result


def get_f_vector(functions_list, first_solution):
    dimension = len(functions_list)
    result = [0 for _ in range(dimension)]

    for i in range(dimension):
        result[i] = functions_list[i](first_solution)

    return result


def vector_norm(vector):
    dimension = len(vector)
    result = 0

    for i in range(dimension):
        result += vector[i] * vector[i]

    return result**(0.5)


def inverse_interpolation_helper(x1, x2, x3, y1, y2, y3):
    a = (y2*y3*x1) / ((y1 - y2)*(y1 - y3))
    b = (y1*y3*x2) / ((y2 - y1)*(y2 - y3))
    c = (y1*y2*x3) / ((y3 - y1)*(y3 - y2))
    return a + b + c


def derivate(function, value):
    delta = 10**(-10)
    numerator = function(value + delta) - function(value)
    denominator = delta
    result = numerator/denominator

    return result


def get_inverse_matrix(matrix):
    cofactors = []
    determinant = matrix_determinant(matrix)
    if(determinant == 0):
        return 0

    for r in range(len(matrix)):
        cofactorRow = []

        for c in range(len(matrix)):
            minor = inverse_auxiliar_function(matrix, r, c)
            cofactorRow.append(((-1)**(r+c)) * matrix_determinant(minor))

        cofactors.append(cofactorRow)

    cofactors = get_transposed_matrix(cofactors)

    for r in range(len(cofactors)):
        for c in range(len(cofactors)):

            cofactors[r][c] = cofactors[r][c]/determinant

    return cofactors


def inverse_auxiliar_function(matrix, i, j):
    return [row[:j] + row[j+1:] for row in (matrix[:i]+matrix[i+1:])]


def sum_vectors(vector_a, vector_b):
    dimension = len(vector_a)
    result = [0 for i in range(dimension)]

    for i in range(dimension):
        result[i] = vector_a[i] + vector_b[i]

    return result


def sum_matrixes(matriz_a, matriz_b):
    dimension = len(matriz_a)
    result = [[0 for i in range(dimension)] for j in range(dimension)]

    for i in range(dimension):
        for j in range(dimension):
            result[i][j] = matriz_a[i][j] + matriz_b[i][j]

    return result


def broyden_method_helper(vector_a, vector_b):
    dimension = len(vector_a)
    result = [[0 for i in range(dimension)] for j in range(dimension)]

    for i in range(dimension):
        for j in range(dimension):
            result[i][j] = vector_a[i]*vector_b[j]

    return result


def get_biggest_element(matrix_a):
    number_of_rows = len(matrix_a)
    biggest_value = -math.inf
    for i in range(number_of_rows):
        for j in range(number_of_rows):
            if (i != j and math.fabs(matrix_a[i][j]) > biggest_value):
                biggest_value = math.fabs(matrix_a[i][j])
                index = (i, j)
    return index


def multiply_matrix_vector(matrix_a, vector):
    number_of_columns_of_a = len(matrix_a)
    number_of_columns_of_vector = range(len(vector))
    result = [0.0 for _ in range(number_of_columns_of_a)]

    for j in range(number_of_columns_of_a):
        summation = 0

        for i in number_of_columns_of_vector:
            summation += matrix_a[j][i]*vector[i]

        result[j] = summation

    return result


def multiply_matrix_scalar(matrix_a, scalar):
    n = len(matrix_a)
    result = [[0.0 for _ in range(n)] for _ in range(n)]

    for j in range(n):
        for i in range(n):
            result[i][j] = matrix_a[i][j] * scalar

    return result


def multiply_matrixes(matrix_a, matrix_b):
    number_of_rows = len(matrix_a)
    result = [[0.0 for _ in range(number_of_rows)]
              for _ in range(number_of_rows)]
    for i in range(len(matrix_a)):
        for j in range(len(matrix_b[0])):
            for k in range(len(matrix_b)):
                result[i][j] += matrix_a[i][k] * matrix_b[k][j]
    return result


def display_matrix(matrix):
    print('\n'.join([''.join(['{:4}'.format(item)
          for item in row]) for row in matrix]))

    print('\n')


def matrix_determinant(matrix):
    number_of_rows = len(matrix)
    number_of_columns = len(matrix[0])

    result = 0

    if(number_of_rows != number_of_columns):
        return "Não é possível calcular a determinante de uma matriz não quadrada"

    if(number_of_rows == 1):
        return matrix[0][0]

    for k in range(number_of_rows):
        result += matrix[k][0]*((-1)**k) * \
            matrix_determinant(get_auxiliar_matrix(matrix, k))

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
    number_of_rows = len(matrix)
    number_of_columns = len(matrix[0])

    if(number_of_rows != number_of_columns):
        return -1

    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if(matrix[i][j] != matrix[j][i]):
                return False

    return True


def sylvesters_criterion(matrix):
    for i in range(len(matrix)):
        secondary = get_main_minor(matrix, i)
        if(matrix_determinant(secondary) <= 0):
            return False

    return True


def get_main_minor(matrix, index):
    result = [[matrix[row][column]
               for row in range(index + 1)] for column in range(index + 1)]
    return result


def backward_substitution(matrix_u, matrix_y):
    number_of_rows = len(matrix_u)
    matrix_x = [0 for i in range(number_of_rows)]

    matrix_x[number_of_rows-1] = matrix_y[number_of_rows-1] / \
        matrix_u[number_of_rows-1][number_of_rows-1]

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
        lines_summation = 0
        columns_summation = 0
        for j in range(len(matrix)):
            if (i != j):
                lines_summation += math.fabs(matrix[i][j])
                columns_summation += math.fabs(matrix[j][i])

        if(matrix[i][i] < lines_summation or matrix[i][i] < columns_summation):
            return False

    return True


def vector_multiplication(vector_a, vector_b):
    result = 0

    for i in range(len(vector_a)):
        for j in range(len(vector_b)):
            if (i == j):
                result += vector_a[i]*vector_b[i]

    return result


def get_transposed_matrix(matrix):
    answer = [[0.0]*len(matrix) for i in range(len(matrix[0]))]

    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            answer[j][i] = matrix[i][j]

    return answer


def phi(matrix, indexes):
    denominator = (matrix[indexes[0]][indexes[0]] -
                   matrix[indexes[1]][indexes[1]])

    if(matrix[indexes[0]][indexes[0]] == matrix[indexes[1]][indexes[1]]):
        return math.pi/4
    else:
        return math.atan(2*matrix[indexes[0]][indexes[1]]/denominator)/2


def calculate_p_matrix(matrix, indexes):
    number_of_rows = len(matrix)
    phi_value = 0

    p_matrix = [[0.0 for _ in range(number_of_rows)]
                for _ in range(number_of_rows)]
    for i in range(number_of_rows):
        p_matrix[i][i] = 1.0

    phi_value = phi(matrix, indexes)

    p_matrix[indexes[0]][indexes[0]] = math.cos(phi_value)
    p_matrix[indexes[1]][indexes[1]] = math.cos(phi_value)
    p_matrix[indexes[0]][indexes[1]] = -math.sin(phi_value)
    p_matrix[indexes[1]][indexes[0]] = math.sin(phi_value)

    return p_matrix


def build_mmq_p_matrix(vector_x):
    number_of_rows = len(vector_x)
    result_matrix = [[1 for x in range(2)] for y in range(number_of_rows)]
    for i in range(number_of_rows):
        result_matrix[i][1] = vector_x[i]

    return result_matrix
