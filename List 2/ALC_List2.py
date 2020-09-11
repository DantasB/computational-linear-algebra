import Matrix_Utils
import math


def power_method(matrix):
    number_of_rows = len(matrix)
    eigenvector    = [1.0 for _ in range(number_of_rows)]

    y_vector      = Matrix_Utils.multiply_matrix_vector(matrix, eigenvector)
    first_element = 1
    eigenvalue    = y_vector[0]
    steps             = 1

    for i in range(number_of_rows):
        y_vector[i] = y_vector[i]/eigenvalue
        eigenvector = y_vector

    residue = math.fabs(eigenvalue - first_element)/eigenvalue
    tol     = 10**-3

    while (residue>=tol):
        first_element = eigenvalue
        y_vector      = Matrix_Utils.multiply_matrix_vector(matrix, eigenvector)
        eigenvalue    = y_vector[0]

        for i in range(number_of_rows):
            y_vector[i] = y_vector[i]/eigenvalue
            eigenvector = y_vector

        residue = math.fabs(eigenvalue-first_element)/eigenvalue

        steps+=1

    print("Eigenvalue: " + str(eigenvalue))
    print("Eigenvector: " + str(eigenvector))
    print("Steps: " + str(steps))


def jacobi_method(matrix):
    if(not Matrix_Utils.check_simetry(matrix)):
        return -1

    number_of_rows  = len(matrix)
    identity_matrix = [[float(i==j) for j in range(number_of_rows)] for i in range(number_of_rows)]
    tol             = 10**(-2)
    biggest_element = Matrix_Utils.get_biggest_element(matrix)

    while (matrix[biggest_element[0]][biggest_element[1]] > tol):
        p_matrix            = Matrix_Utils.calculate_p_matrix(matrix, biggest_element)
        p_matrix_transposed = Matrix_Utils.get_transposed_matrix(p_matrix)
        matrix              = Matrix_Utils.multiply_matrixes(p_matrix_transposed, Matrix_Utils.multiply_matrixes(matrix, p_matrix))
        identity_matrix     = Matrix_Utils.multiply_matrixes(identity_matrix, p_matrix)
        biggest_element     = Matrix_Utils.get_biggest_element(matrix)

    autovalores = []
    for i in range(number_of_rows):
        autovalores.append(matrix[i][i])

    print("Autovalores: " + str(autovalores))
    print("Autovetores: " + str(identity_matrix))
