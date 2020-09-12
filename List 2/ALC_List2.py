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
    tol             = 10**(-3)
    biggest_element = Matrix_Utils.get_biggest_element(matrix)
    steps           = 1

    while (math.fabs(matrix[biggest_element[0]][biggest_element[1]]) > tol):
        p_matrix            = Matrix_Utils.calculate_p_matrix(matrix, biggest_element)
        p_matrix_transposed = Matrix_Utils.get_transposed_matrix(p_matrix)
        matrix              = Matrix_Utils.multiply_matrixes(p_matrix_transposed, Matrix_Utils.multiply_matrixes(matrix, p_matrix))
        identity_matrix     = Matrix_Utils.multiply_matrixes(identity_matrix, p_matrix)
        biggest_element     = Matrix_Utils.get_biggest_element(matrix)

        steps+=1
    result = []
    for i in range(number_of_rows):
            result.append((i+1, matrix[i][i], identity_matrix[i]))
            print(str(i+1) + "° Autovalor: " + str(matrix[i][i]) + ", Autovetor: " + str(identity_matrix[i]))
            
    print("Steps: "+ str(steps))
    
    return result
