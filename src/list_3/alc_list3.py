from src.utils.matrix_utils import build_mmq_p_matrix, get_transposed_matrix, multiply_matrixes, multiply_matrix_vector
from src.list_1.alc_list1 import solve_matrix


def mmq(vector_x, vector_y):
    matrix_p = build_mmq_p_matrix(vector_x)
    matrix_p_transposed = get_transposed_matrix(matrix_p)
    matrix_a = multiply_matrixes(matrix_p_transposed, matrix_p)
    matrix_c = multiply_matrix_vector(
        matrix_p_transposed, vector_y)
    return solve_matrix(matrix_a, matrix_c, False)
