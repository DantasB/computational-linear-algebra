from importlib.machinery import SourceFileLoader

Matrix_Utils = SourceFileLoader("Matrix_Utils", "/home/bdantas/Documentos/Repos/Algebra-Linear-Computacional/Utils/Matrix_Utils.py").load_module()
ALC_List1    = SourceFileLoader("ALC_List1", "/home/bdantas/Documentos/Repos/Algebra-Linear-Computacional/List 1/ALC_List1.py").load_module()

def mmq(vector_x, vector_y):
    matrix_p            = Matrix_Utils.build_mmq_p_matrix(vector_x);
    matrix_p_transposed = Matrix_Utils.get_transposed_matrix(matrix_p)
    matrix_a            = Matrix_Utils.multiply_matrixes(matrix_p_transposed, matrix_p)
    matrix_c            = Matrix_Utils.multiply_matrix_vector(matrix_p_transposed, vector_y)
    return ALC_List1.solve_matrix(matrix_a, matrix_c, False);
