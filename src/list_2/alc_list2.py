import math

from src.utils.matrix_utils import multiply_matrix_vector, check_simetry, get_biggest_element, calculate_p_matrix, get_transposed_matrix, multiply_matrixes


def power_method(matrix):
    number_of_rows = len(matrix)
    eigenvector = [1.0 for _ in range(number_of_rows)]

    y_vector = multiply_matrix_vector(matrix, eigenvector)
    first_element = 1
    eigenvalue = y_vector[0]
    steps = 1

    for i in range(number_of_rows):
        y_vector[i] = y_vector[i]/eigenvalue
        eigenvector = y_vector

    residue = math.fabs(eigenvalue - first_element)/eigenvalue
    tol = 10**-5

    while (residue >= tol):
        first_element = eigenvalue
        y_vector = multiply_matrix_vector(matrix, eigenvector)
        eigenvalue = y_vector[0]

        for i in range(number_of_rows):
            y_vector[i] = y_vector[i]/eigenvalue

        eigenvector = y_vector
        residue = math.fabs(eigenvalue-first_element)/eigenvalue

        steps += 1

    print("Eigenvalue: " + str(eigenvalue))
    print("Eigenvector: " + str(eigenvector))
    print("Steps: " + str(steps))


def jacobi_method(matrix):
    if(not check_simetry(matrix)):
        return -1

    number_of_rows = len(matrix)
    eigenvector = [[float(i == j) for j in range(number_of_rows)]
                   for i in range(number_of_rows)]
    tol = 10**(-5)
    biggest_element = get_biggest_element(matrix)
    steps = 1

    while (math.fabs(matrix[biggest_element[0]][biggest_element[1]]) > tol):
        p_matrix = calculate_p_matrix(matrix, biggest_element)
        p_matrix_transposed = get_transposed_matrix(p_matrix)
        matrix = multiply_matrixes(
            p_matrix_transposed, multiply_matrixes(matrix, p_matrix))
        eigenvector = multiply_matrixes(eigenvector, p_matrix)
        biggest_element = get_biggest_element(matrix)

        steps += 1
    result = []
    for i in range(number_of_rows):
        result.append((i+1, matrix[i][i], eigenvector[i]))
        print(str(i+1) + "° Eigenvalue: " +
              str(matrix[i][i]) + ", Eigenvector: " + str(eigenvector[i]))

    print("Steps: " + str(steps))

    return result
